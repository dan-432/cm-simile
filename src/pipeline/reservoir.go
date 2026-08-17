package pipeline

import (
	"context"
	"fmt"
	"math"
	"math/rand"
	"sort"

	"cmsimile/src/matrix"
	"cmsimile/src/seqpool"
	"cmsimile/src/shuffle"
)

// ReservoirConfig controls the adaptive null-reservoir construction. See
// the package doc comment and buildReservoir's doc comment for the full
// design: a reservoir is a large, once-scored pool of dinucleotide-
// shuffled sequences per model, grown until a resampling diagnostic
// stabilizes, then reused via bootstrap resampling for every pair that
// model participates in.
type ReservoirConfig struct {
	// ResampleSize is the size of each bootstrap draw used both for
	// convergence diagnostics and for the final per-pair null replicates.
	// Defaults to NPool if zero -- keeping null replicate means computed
	// over the same sample size as the real pool's own mean, so the two
	// are on a comparable footing.
	ResampleSize int

	// ReservoirMin is the starting reservoir size per model, and the
	// floor below which ResampleSize must not exceed it (you can't
	// resample-with-replacement usefully from fewer distinct values than
	// you're drawing, in the sense that variance estimates below this
	// point are unreliable -- see validate()).
	ReservoirMin int

	// ReservoirMax caps how large a reservoir is allowed to grow. If a
	// model hits this cap without its diagnostic converging (within
	// ReservoirTolerance), growth stops anyway and a warning is surfaced
	// via Warnings on the result -- silently proceeding with an
	// unconverged reservoir would be exactly the kind of quiet
	// miscalibration this whole redesign is trying to avoid repeating.
	ReservoirMax int

	// ReservoirTolerance is the maximum relative change in the
	// convergence diagnostic (see diagnosticTailQuantile) between one
	// doubling round and the next that counts as "converged".
	ReservoirTolerance float64

	// DiagnosticResamples is how many bootstrap resamples are drawn to
	// *estimate* the convergence diagnostic each round -- deliberately a
	// separate, usually smaller, number than NShuffles (the final
	// per-pair null replicate count), since the diagnostic only needs to
	// be stable enough to detect convergence, not to be the final
	// reported null distribution itself.
	DiagnosticResamples int
}

func (rc ReservoirConfig) resampleSize(npool int) int {
	if rc.ResampleSize > 0 {
		return rc.ResampleSize
	}
	return npool
}

// reservoirResult holds the fully-grown reservoir (pool + score matrix)
// plus bookkeeping surfaced back to the caller.
type reservoirResult struct {
	Pool     *seqpool.Pool
	Matrix   *matrix.ScoreMatrix
	Sizes    map[string]int // final reservoir size per model
	Rounds   map[string]int // number of growth rounds each model needed
	Warnings []string       // e.g. "model X hit ReservoirMax without converging"
}

// buildReservoir grows one shared, multi-origin reservoir of
// dinucleotide-shuffled sequences -- doubling in size each round,
// starting at rc.ReservoirMin -- until every model's convergence
// diagnostic (see diagnosticTailQuantile) stabilizes within
// rc.ReservoirTolerance, or rc.ReservoirMax is hit.
//
// Each round scores only the newly-added increment (not the whole
// reservoir again) against every model in cfg.Models, via the same
// scorePoolAgainstAllModels batching used for the real pool -- so total
// cmsearch calls are O(NModels x log(finalSize/ReservoirMin)), not
// O(NModels x NShuffles) the way the old per-replicate shuffle-and-
// rescore design was. This is what makes NShuffles cheap afterward: it
// now only controls how many bootstrap resamples are drawn from an
// already-scored reservoir, in memory, with zero further subprocess
// calls.
//
// A model only needs to converge against the OTHER models it actually
// gets compared to (not against itself -- the self-baseline in
// stats.ComputePairStat always comes from the real pool, never the
// reservoir, so a model's own reservoir is never scored-and-used against
// itself). A model that converges early stops contributing new sequences
// to later rounds' increment pool, but its already-scored reservoir
// entries remain part of the final reservoir.
func buildReservoir(ctx context.Context, cfg Config, realPool *seqpool.Pool, rc ReservoirConfig) (*reservoirResult, error) {
	if err := rc.validate(cfg.NPool); err != nil {
		return nil, err
	}

	reservoirPool := &seqpool.Pool{}
	reservoirMatrix := matrix.New()

	size := make(map[string]int, len(cfg.Models))
	converged := make(map[string]bool, len(cfg.Models))
	rounds := make(map[string]int, len(cfg.Models))
	prevDiag := make(map[string]map[string]float64, len(cfg.Models)) // model -> target -> diagnostic
	var warnings []string

	targetSize := rc.ReservoirMin
	round := 0
	for {
		growing := make([]Model, 0, len(cfg.Models))
		for _, m := range cfg.Models {
			if !converged[m.Name] {
				growing = append(growing, m)
			}
		}
		if len(growing) == 0 {
			break
		}

		increment := &seqpool.Pool{}
		for _, m := range growing {
			need := targetSize - size[m.Name]
			if need <= 0 {
				continue
			}
			templates := realPool.ByOrigin(m.Name)
			if len(templates) == 0 {
				return nil, fmt.Errorf("pipeline: no real pool sequences to shuffle for model %q", m.Name)
			}
			rng := rand.New(rand.NewSource(cfg.Seed + reservoirSeedOffset(m.Name, round)))
			for k := 0; k < need; k++ {
				template := templates[k%len(templates)]
				shuffled, err := shuffle.DinucleotideShuffle(template.Residues, rng)
				if err != nil {
					return nil, fmt.Errorf("pipeline: shuffling reservoir sequence for %q: %w", m.Name, err)
				}
				id := fmt.Sprintf("%s_reservoir_%06d", m.Name, size[m.Name]+k)
				increment.Add(seqpool.Sequence{ID: id, Origin: m.Name, Residues: shuffled})
			}
			size[m.Name] = targetSize
			rounds[m.Name]++
		}

		if increment.Len() == 0 {
			break
		}

		stage := fmt.Sprintf("scoring_reservoir_round_%d", round)
		tick := newTicker(cfg.OnProgress, stage, len(cfg.Models))
		incrementMatrix, err := scorePoolAgainstAllModels(ctx, cfg, increment, tick)
		if err != nil {
			return nil, fmt.Errorf("pipeline: scoring reservoir growth round %d: %w", round, err)
		}

		reservoirPool.Add(increment.Sequences...)
		for _, seq := range increment.Sequences {
			for _, m := range cfg.Models {
				sc, ok := incrementMatrix.Get(seq.ID, m.Name)
				if ok {
					reservoirMatrix.Set(seq.ID, m.Name, sc)
				}
			}
		}

		// Convergence check, per growing model, against every other
		// model it will actually be compared to. A model only converges
		// once it's had at least two rounds' worth of diagnostic values
		// to compare (prevDiag existing is what signals "not the first
		// measurement").
		diagRNG := rand.New(rand.NewSource(cfg.Seed + 3_000_017 + int64(round)*13))
		for _, m := range growing {
			ids := idsOf(reservoirPool.ByOrigin(m.Name))
			diagNow := make(map[string]float64, len(cfg.Models)-1)
			for _, target := range cfg.Models {
				if target.Name == m.Name {
					continue // reservoir is never scored-and-used against its own origin -- see doc comment
				}
				diagNow[target.Name] = diagnosticTailQuantile(reservoirMatrix, ids, target.Name, rc.resampleSize(cfg.NPool), rc.DiagnosticResamples, diagRNG)
			}

			if prev, hadPrev := prevDiag[m.Name]; hadPrev {
				worst := 0.0
				for target, now := range diagNow {
					worst = math.Max(worst, relativeChange(now, prev[target]))
				}
				if worst < rc.ReservoirTolerance {
					converged[m.Name] = true
				}
			}
			prevDiag[m.Name] = diagNow

			if !converged[m.Name] && size[m.Name] >= rc.ReservoirMax {
				converged[m.Name] = true
				warnings = append(warnings, fmt.Sprintf(
					"model %q: reservoir hit ReservoirMax (%d) after %d rounds without the convergence diagnostic stabilizing within tolerance %.4f -- null p-values for pairs involving this model rest on a less-validated reservoir size than the rest of the run",
					m.Name, rc.ReservoirMax, rounds[m.Name], rc.ReservoirTolerance))
			}
		}

		targetSize *= 2
		round++
	}

	return &reservoirResult{
		Pool:     reservoirPool,
		Matrix:   reservoirMatrix,
		Sizes:    size,
		Rounds:   rounds,
		Warnings: warnings,
	}, nil
}

func (rc ReservoirConfig) validate(npool int) error {
	if rc.ReservoirMin <= 0 {
		return fmt.Errorf("pipeline: ReservoirConfig.ReservoirMin must be > 0")
	}
	if rc.ReservoirMax < rc.ReservoirMin {
		return fmt.Errorf("pipeline: ReservoirConfig.ReservoirMax (%d) must be >= ReservoirMin (%d)", rc.ReservoirMax, rc.ReservoirMin)
	}
	if rc.ReservoirTolerance <= 0 {
		return fmt.Errorf("pipeline: ReservoirConfig.ReservoirTolerance must be > 0")
	}
	if rc.DiagnosticResamples <= 0 {
		return fmt.Errorf("pipeline: ReservoirConfig.DiagnosticResamples must be > 0")
	}
	resampleSize := rc.resampleSize(npool)
	if resampleSize > rc.ReservoirMin {
		return fmt.Errorf("pipeline: resample size (%d) must not exceed ReservoirMin (%d) -- "+
			"can't draw a meaningful bootstrap sample larger than the starting reservoir", resampleSize, rc.ReservoirMin)
	}
	return nil
}

// diagnosticTailQuantile bootstrap-resamples resamples draws of size
// resampleSize (with replacement) from ids' scores against target, and
// returns the 5th percentile of the resulting means. The tail, not the
// center, is what a p-value threshold actually depends on, so that's
// what's tracked for convergence rather than e.g. the plain mean or
// variance -- a reservoir whose central tendency has stabilized can
// still have an unstable tail if it's too small.
func diagnosticTailQuantile(m *matrix.ScoreMatrix, ids []string, target string, resampleSize, nResamples int, rng *rand.Rand) float64 {
	medians := bootstrapMedians(m, ids, target, resampleSize, nResamples, rng)
	sort.Float64s(medians)
	idx := int(0.05 * float64(len(medians)))
	if idx >= len(medians) {
		idx = len(medians) - 1
	}
	return medians[idx]
}

// bootstrapMeans draws nResamples bootstrap samples (with replacement) of
// size resampleSize from ids' scores against target, returning each
// sample's mean. This is the core in-memory resampling primitive used
// both for the convergence diagnostic above and for the final per-pair
// null replicates in pipeline.go -- no cmsearch calls, pure arithmetic
// over an already-scored reservoir.
func bootstrapMeans(m *matrix.ScoreMatrix, ids []string, target string, resampleSize, nResamples int, rng *rand.Rand) []float64 {
	means := make([]float64, nResamples)
	for i := 0; i < nResamples; i++ {
		var sum float64
		for j := 0; j < resampleSize; j++ {
			id := ids[rng.Intn(len(ids))]
			sc, _ := m.Get(id, target)
			sum += sc
		}
		means[i] = sum / float64(resampleSize)
	}
	return means
}

// bootstrapMedians draws nResamples bootstrap samples (with replacement) of
// size resampleSize from ids' scores against target, returning each
// sample's median. This is the core in-memory resampling primitive used
// both for the convergence diagnostic above and for the final per-pair
// null replicates in pipeline.go -- no cmsearch calls, pure arithmetic
// over an already-scored reservoir.
func bootstrapMedians(m *matrix.ScoreMatrix, ids []string, target string, resampleSize, nResamples int, rng *rand.Rand) []float64 {
	medians := make([]float64, nResamples)
	for i := 0; i < nResamples; i++ {
		sample := make([]float64, resampleSize)
		for j := 0; j < resampleSize; j++ {
			id := ids[rng.Intn(len(ids))]
			sc, _ := m.Get(id, target)
			sample[j] = sc
		}
		sort.Float64s(sample)
		medians[i] = sample[resampleSize/2]
	}
	return medians
}

func relativeChange(now, prev float64) float64 {
	denom := math.Abs(prev)
	if denom < 1e-9 {
		denom = 1e-9
	}
	return math.Abs(now-prev) / denom
}

// reservoirSeedOffset gives each (model, round) combination its own
// deterministic sub-seed, distinct from the offsets used elsewhere
// (buildRealPool, the old per-replicate shuffle seeding) so reservoir
// growth doesn't accidentally reuse a random stream from another stage.
func reservoirSeedOffset(modelName string, round int) int64 {
	var h int64 = 5_381
	for _, c := range modelName {
		h = h*33 + int64(c)
	}
	return 9_000_011 + h%1_000_000 + int64(round)*1_000_003
}

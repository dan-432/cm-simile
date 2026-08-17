// Package pipeline orchestrates the full cm-simile run: build one real
// pool across all input models, score it once per model, build an
// adaptively-sized reservoir of dinucleotide-shuffled null sequences and
// score that once too, then generate every pair's null distribution via
// cheap in-memory bootstrap resampling from the reservoir (see
// reservoir.go) and slice all pairwise statistics out of the resulting
// matrices.
//
// This is deliberately N-way from the ground up rather than pairwise-first:
// 1-to-1, bipartite (query set vs. reference set), and full many-to-many
// comparison are all the same Run() call, differing only in how many
// Models are supplied and which pairs Query restricts the *output* to.
// Scoring cost still grows with the number of models actually supplied
// (see the README's "Algorithm Approach" section on the bipartite vs.
// full-all-vs-all cost distinction) -- Query filters which pairs are
// reported, not which models get pooled and scored.
package pipeline

import (
	"context"
	"fmt"
	"math/rand"
	"sort"
	"sync"
	"sync/atomic"

	"cmsimile/src/matrix"
	"cmsimile/src/seqpool"
	"cmsimile/src/stats"
)

// Model identifies one covariance model to include in the comparison.
type Model struct {
	Name string // used as the model's identifier throughout output
	Path string // path to the .cm file
}

// Emitter is the subset of infernal.Emitter's method set the pipeline
// needs. Defined here (rather than depending on the infernal package
// directly) so tests can supply a fake without needing Infernal installed;
// infernal.Emitter satisfies this interface structurally.
type Emitter interface {
	Sample(ctx context.Context, cmPath, originLabel string, n int, seed int64) ([]seqpool.Sequence, error)
}

// Scorer is the subset of infernal.Scorer's method set the pipeline needs.
// infernal.Scorer satisfies this interface structurally.
type Scorer interface {
	Score(ctx context.Context, cmPath string, seqs []seqpool.Sequence) (map[string]float64, error)
}

// Config controls a pipeline run.
type Config struct {
	Models []Model

	// Query, if non-empty, restricts the *reported* pairs to those where
	// at least one model's name is in Query. All supplied Models are
	// still pooled and scored regardless (see package doc) -- this only
	// filters output, e.g. for a bipartite "a few novel motifs against
	// all of Rfam" run where you only care about pairs touching your
	// query motifs. Leave nil/empty to report every pair.
	Query map[string]bool

	NPool int // real sequences sampled per model

	// NShuffles is now the number of bootstrap resamples drawn per pair
	// from the (once-built, once-scored) null reservoir -- see
	// reservoir.go. This is cheap: resampling is in-memory arithmetic,
	// not a subprocess call, so NShuffles no longer drives cmsearch call
	// count the way it did under the old per-replicate shuffle-and-
	// rescore design. Reservoir construction cost is instead governed by
	// Reservoir's own fields.
	NShuffles int

	Reservoir ReservoirConfig

	Seed    int64 // base seed; each model/replicate derives its own sub-seed from this for reproducibility
	Threads int   // max concurrent subprocess calls

	Emitter Emitter
	Scorer  Scorer

	// OnProgress, if set, is called as work completes: once per cmemit
	// sampling call, once per cmsearch scoring call for the real pool,
	// and once per cmsearch scoring call for each reservoir growth round
	// (stage name includes the round number, e.g.
	// "scoring_reservoir_round_0", since the number of rounds is
	// data-dependent and not known upfront). completed/total are
	// cumulative within that stage. Safe to call from multiple
	// goroutines concurrently -- callers doing anything beyond simple,
	// fast, thread-safe reporting (e.g. printing) should synchronize
	// internally.
	OnProgress func(stage string, completed, total int)

	// DumpMatricesDir, if non-empty, writes the real pool's and the null
	// reservoir's full (sequence x model) score matrices to this
	// directory as long-format TSVs -- e.g. for plotting the null
	// distribution's shape directly, rather than only its summary
	// statistics. See DumpMatrices in matrixdump.go.
	DumpMatricesDir string
}

// WorkEstimate summarises how much subprocess-level work a Config implies,
// for printing an overview before a run starts. See the package doc
// comment and the top-level README's complexity discussion for why
// subprocess *call count* (linear in model count) and total *scoring
// compute* (quadratic, since each call's pool grows with model count) are
// different things worth showing separately.
//
// Reservoir cost is a RANGE, not an exact number, because reservoir.go's
// growth is adaptive -- it stops once a convergence diagnostic
// stabilizes, which depends on the actual score data, not just the
// config. CmsearchCallsReservoirBest assumes every model converges after
// the first round (ReservoirMin sequences, no growth needed);
// CmsearchCallsReservoirWorst assumes every model grows all the way to
// ReservoirMax without converging. Real runs should land somewhere
// between the two -- printing a false-precision single number here would
// misrepresent what's actually a data-dependent quantity.
type WorkEstimate struct {
	NModels        int
	NPairsTotal    int // all pairs among NModels
	NPairsReported int // pairs surviving the Query filter
	CmemitCalls    int // one per model

	CmsearchCallsReal int // one per model, scoring the real pool

	CmsearchCallsReservoirBest  int // best case: every model converges in round 1
	CmsearchCallsReservoirWorst int // worst case: every model grows to ReservoirMax

	CmsearchCallsTotalBest  int
	CmsearchCallsTotalWorst int
}

// Estimate computes a WorkEstimate for cfg without running anything.
func Estimate(cfg Config) WorkEstimate {
	n := len(cfg.Models)
	totalPairs := n * (n - 1) / 2
	reported := 0
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			if len(cfg.Query) > 0 && !cfg.Query[cfg.Models[i].Name] && !cfg.Query[cfg.Models[j].Name] {
				continue
			}
			reported++
		}
	}
	cmsearchReal := n

	best := n // one round, ReservoirMin sequences, every model converges immediately
	worst := 0
	if cfg.Reservoir.ReservoirMin > 0 && cfg.Reservoir.ReservoirMax >= cfg.Reservoir.ReservoirMin {
		rounds := 1
		size := cfg.Reservoir.ReservoirMin
		for size < cfg.Reservoir.ReservoirMax {
			size *= 2
			rounds++
		}
		worst = n * rounds
	}

	return WorkEstimate{
		NModels:                     n,
		NPairsTotal:                 totalPairs,
		NPairsReported:              reported,
		CmemitCalls:                 n,
		CmsearchCallsReal:           cmsearchReal,
		CmsearchCallsReservoirBest:  best,
		CmsearchCallsReservoirWorst: worst,
		CmsearchCallsTotalBest:      cmsearchReal + best,
		CmsearchCallsTotalWorst:     cmsearchReal + worst,
	}
}

// Run executes the full pipeline and returns one stats.PairStat per
// reported pair (sorted by ascending min(p-value), most significant
// overlap first, with model-name pairs as a stable tie-break), plus any
// warnings surfaced during reservoir construction (see
// reservoirResult.Warnings) that are worth a human looking at even
// though they didn't stop the run.
func Run(ctx context.Context, cfg Config) ([]stats.PairStat, []string, error) {
	if err := validate(cfg); err != nil {
		return nil, nil, err
	}

	realPool, err := buildRealPool(ctx, cfg)
	if err != nil {
		return nil, nil, fmt.Errorf("pipeline: building real pool: %w", err)
	}

	realMatrix, err := scorePoolAgainstAllModels(ctx, cfg, realPool, newTicker(cfg.OnProgress, "scoring_real", len(cfg.Models)))
	if err != nil {
		return nil, nil, fmt.Errorf("pipeline: scoring real pool: %w", err)
	}

	reservoir, err := buildReservoir(ctx, cfg, realPool, cfg.Reservoir)
	if err != nil {
		return nil, nil, fmt.Errorf("pipeline: building null reservoir: %w", err)
	}

	if cfg.DumpMatricesDir != "" {
		if err := DumpMatrices(cfg.DumpMatricesDir, realPool, realMatrix, reservoir.Pool, reservoir.Matrix); err != nil {
			return nil, nil, fmt.Errorf("pipeline: dumping matrices: %w", err)
		}
	}

	results, err := computeAllPairStats(cfg, realPool, realMatrix, reservoir)
	if err != nil {
		return nil, nil, fmt.Errorf("pipeline: computing pairwise statistics: %w", err)
	}

	sort.Slice(results, func(i, j int) bool {
		pi := minP(results[i])
		pj := minP(results[j])
		if pi != pj {
			return pi < pj
		}
		if results[i].ModelA != results[j].ModelA {
			return results[i].ModelA < results[j].ModelA
		}
		return results[i].ModelB < results[j].ModelB
	})
	return results, reservoir.Warnings, nil
}

func validate(cfg Config) error {
	if len(cfg.Models) < 2 {
		return fmt.Errorf("pipeline: need at least 2 models, got %d", len(cfg.Models))
	}
	seen := make(map[string]bool, len(cfg.Models))
	for _, m := range cfg.Models {
		if seen[m.Name] {
			return fmt.Errorf("pipeline: duplicate model name %q", m.Name)
		}
		seen[m.Name] = true
	}
	if cfg.NPool <= 0 {
		return fmt.Errorf("pipeline: NPool must be > 0, got %d", cfg.NPool)
	}
	if cfg.NShuffles <= 0 {
		return fmt.Errorf("pipeline: NShuffles must be > 0, got %d", cfg.NShuffles)
	}
	if cfg.Threads <= 0 {
		return fmt.Errorf("pipeline: Threads must be > 0, got %d", cfg.Threads)
	}
	if cfg.Emitter == nil || cfg.Scorer == nil {
		return fmt.Errorf("pipeline: Emitter and Scorer must be set")
	}
	if err := cfg.Reservoir.validate(cfg.NPool); err != nil {
		return err
	}
	return nil
}

// progressTicker is a small thread-safe counter used to drive
// Config.OnProgress from concurrent goroutines without each call site
// needing its own synchronization.
type progressTicker struct {
	stage string
	total int
	done  int64
	fn    func(stage string, completed, total int)
}

func newTicker(fn func(stage string, completed, total int), stage string, total int) *progressTicker {
	if fn == nil {
		return nil
	}
	return &progressTicker{stage: stage, total: total, fn: fn}
}

func (t *progressTicker) tick() {
	if t == nil {
		return
	}
	n := atomic.AddInt64(&t.done, 1)
	t.fn(t.stage, int(n), t.total)
}

// buildRealPool draws cfg.NPool sequences from every model, in parallel up
// to cfg.Threads, and combines them into one pool.
func buildRealPool(ctx context.Context, cfg Config) (*seqpool.Pool, error) {
	pool := &seqpool.Pool{}
	var mu sync.Mutex
	tick := newTicker(cfg.OnProgress, "sampling", len(cfg.Models))

	err := parallelFor(len(cfg.Models), cfg.Threads, func(i int) error {
		m := cfg.Models[i]
		// Each model gets its own deterministic sub-seed derived from the
		// base seed and its index, so re-running with the same base seed
		// reproduces the same pools regardless of goroutine scheduling
		// order.
		subSeed := cfg.Seed + int64(i)*1_000_003
		seqs, err := cfg.Emitter.Sample(ctx, m.Path, m.Name, cfg.NPool, subSeed)
		if err != nil {
			return fmt.Errorf("sampling model %q: %w", m.Name, err)
		}
		mu.Lock()
		pool.Add(seqs...)
		mu.Unlock()
		tick.tick()
		return nil
	})
	if err != nil {
		return nil, err
	}
	return pool, nil
}

// scorePoolAgainstAllModels scores the given pool against every model in
// cfg.Models, in parallel up to cfg.Threads, returning the combined
// ScoreMatrix. This is the O(N) subprocess-call step: one Score call per
// model, regardless of how many pairs are ultimately reported. tick, if
// non-nil, is called once per completed model score (used to drive
// Config.OnProgress; shared across shuffle replicates by
// buildNullDistributions so the null-scoring stage reports cumulative
// progress across all replicates, not just the current one).
func scorePoolAgainstAllModels(ctx context.Context, cfg Config, pool *seqpool.Pool, tick *progressTicker) (*matrix.ScoreMatrix, error) {
	m := matrix.New()
	var mu sync.Mutex

	err := parallelFor(len(cfg.Models), cfg.Threads, func(i int) error {
		model := cfg.Models[i]
		scores, err := cfg.Scorer.Score(ctx, model.Path, pool.Sequences)
		if err != nil {
			return fmt.Errorf("scoring pool against model %q: %w", model.Name, err)
		}
		mu.Lock()
		m.SetAll(model.Name, scores)
		mu.Unlock()
		tick.tick()
		return nil
	})
	if err != nil {
		return nil, err
	}
	return m, nil
}

// (buildNullDistributions and shuffleWholePool, which used to live here,
// are replaced entirely by reservoir.go's buildReservoir -- see the
// package doc comment.)

// computeAllPairStats slices the real matrix and bootstrap-resamples null
// replicates from the reservoir for each reported pair. Resampling is
// pure in-memory arithmetic (bootstrapMeans in reservoir.go) -- this is
// the payoff of the reservoir redesign: no cmsearch calls happen in this
// function at all, regardless of NShuffles.
func computeAllPairStats(cfg Config, realPool *seqpool.Pool, realMatrix *matrix.ScoreMatrix, reservoir *reservoirResult) ([]stats.PairStat, error) {
	resampleSize := cfg.Reservoir.resampleSize(cfg.NPool)
	var results []stats.PairStat
	for i := 0; i < len(cfg.Models); i++ {
		for j := i + 1; j < len(cfg.Models); j++ {
			a, b := cfg.Models[i].Name, cfg.Models[j].Name
			if len(cfg.Query) > 0 && !cfg.Query[a] && !cfg.Query[b] {
				continue // neither side of this pair is a queried model
			}

			// Deterministic, pair-specific RNG seed so results are
			// reproducible given the same base Seed, independent of any
			// map/goroutine iteration order elsewhere in the pipeline.
			pairSeed := cfg.Seed + 5_000_011 + int64(i)*97 + int64(j)*10_007
			rng := rand.New(rand.NewSource(pairSeed))

			idsA := idsOf(reservoir.Pool.ByOrigin(a))
			idsB := idsOf(reservoir.Pool.ByOrigin(b))
			nullAonB := bootstrapMedians(reservoir.Matrix, idsA, b, resampleSize, cfg.NShuffles, rng)
			nullBonA := bootstrapMedians(reservoir.Matrix, idsB, a, resampleSize, cfg.NShuffles, rng)

			ps, err := stats.ComputePairStat(realPool, realMatrix, a, b, nullAonB, nullBonA)
			if err != nil {
				return nil, fmt.Errorf("pair (%s, %s): %w", a, b, err)
			}
			results = append(results, ps)
		}
	}
	return results, nil
}

func minP(p stats.PairStat) float64 {
	if p.PValueAintoB < p.PValueBintoA {
		return p.PValueAintoB
	}
	return p.PValueBintoA
}

func idsOf(seqs []seqpool.Sequence) []string {
	out := make([]string, len(seqs))
	for i, s := range seqs {
		out[i] = s.ID
	}
	return out
}

// parallelFor runs worker(i) for i in [0, n) using up to `threads`
// concurrent goroutines bounded by a semaphore, and returns the first
// error encountered (in index order, for deterministic error reporting)
// after all launched workers have finished. It does not attempt to
// cancel in-flight subprocess calls on first error -- for this tool's
// scale (comparisons run in minutes, not hours) simplicity and
// correctness matter more than shaving the tail latency of an
// already-failing run.
func parallelFor(n, threads int, worker func(i int) error) error {
	if threads > n {
		threads = n
	}
	if threads < 1 {
		threads = 1
	}
	sem := make(chan struct{}, threads)
	errs := make([]error, n)
	var wg sync.WaitGroup

	for i := 0; i < n; i++ {
		wg.Add(1)
		sem <- struct{}{}
		go func(i int) {
			defer wg.Done()
			defer func() { <-sem }()
			errs[i] = worker(i)
		}(i)
	}
	wg.Wait()

	for _, err := range errs {
		if err != nil {
			return err
		}
	}
	return nil
}

// Package pipeline orchestrates the full cm-simile run: build one real
// pool across all input models, score it once per model, generate
// dinucleotide-shuffled null replicates and score those too, then slice
// all pairwise statistics out of the resulting matrices.
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
	"cmsimile/src/shuffle"
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
	// still pooled and scored regardless - this only
	// filters output, e.g. for a bipartite "a few novel motifs against
	// all of Rfam" run where you only care about pairs touching your
	// query motifs. Leave nil/empty to report every pair.
	Query map[string]bool

	NPool     int   // real sequences sampled per model
	NShuffles int   // dinucleotide-shuffle null replicates
	Seed      int64 // base seed; each model/replicate derives its own sub-seed from this for reproducibility
	Threads   int   // max concurrent subprocess calls

	Emitter Emitter
	Scorer  Scorer

	// OnProgress, if set, is called as work completes: once per cmemit
	// sampling call, and once per cmsearch scoring call (both for the
	// real pool and for each null-shuffle replicate). stage is one of
	// "sampling", "scoring_real", "scoring_null". completed/total are
	// cumulative within that stage across the whole run. Safe to call
	// from multiple goroutines concurrently -- callers doing anything
	// beyond simple, fast, thread-safe reporting (e.g. printing) should
	// synchronize internally.
	OnProgress func(stage string, completed, total int)
}

// WorkEstimate summarises how much subprocess-level work a Config implies,
// for printing an overview before a run starts. See the package doc
// comment and the top-level README's complexity discussion for why
// subprocess *call count* (linear in model count) and total *scoring
// compute* (quadratic, since each call's pool grows with model count) are
// different things worth showing separately.
type WorkEstimate struct {
	NModels            int
	NPairsTotal        int // all pairs among NModels
	NPairsReported     int // pairs surviving the Query filter
	CmemitCalls        int // one per model
	CmsearchCallsReal  int // one per model, scoring the real pool
	CmsearchCallsNull  int // one per model per shuffle replicate
	CmsearchCallsTotal int
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
	cmsearchNull := n * cfg.NShuffles
	return WorkEstimate{
		NModels:            n,
		NPairsTotal:        totalPairs,
		NPairsReported:     reported,
		CmemitCalls:        n,
		CmsearchCallsReal:  cmsearchReal,
		CmsearchCallsNull:  cmsearchNull,
		CmsearchCallsTotal: cmsearchReal + cmsearchNull,
	}
}

// Run executes the full pipeline and returns one stats.PairStat per
// reported pair, sorted by ascending min(p-value) (most significant
// overlap first), with model-name pairs as a stable tie-break.
func Run(ctx context.Context, cfg Config) ([]stats.PairStat, error) {
	if err := validate(cfg); err != nil {
		return nil, err
	}

	realPool, err := buildRealPool(ctx, cfg)
	if err != nil {
		return nil, fmt.Errorf("pipeline: building real pool: %w", err)
	}

	realMatrix, err := scorePoolAgainstAllModels(ctx, cfg, realPool, newTicker(cfg.OnProgress, "scoring_real", len(cfg.Models)))
	if err != nil {
		return nil, fmt.Errorf("pipeline: scoring real pool: %w", err)
	}

	// nullMeanCross[origin][scoredAgainst] = one mean-score value per
	// shuffle replicate: the mean score of origin's shuffled sequences
	// when scored against model `scoredAgainst`.
	nullMeanCross, err := buildNullDistributions(ctx, cfg, realPool)
	if err != nil {
		return nil, fmt.Errorf("pipeline: building null distributions: %w", err)
	}

	results, err := computeAllPairStats(cfg, realPool, realMatrix, nullMeanCross)
	if err != nil {
		return nil, fmt.Errorf("pipeline: computing pairwise statistics: %w", err)
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
	return results, nil
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

// buildNullDistributions generates cfg.NShuffles dinucleotide-shuffled
// replicates of realPool and scores each replicate against every model,
// accumulating per-replicate mean cross-scores.
//
// Subprocess call count: O(NShuffles x len(Models)), the same linear-in-N
// scaling as the real pool, multiplied by the number of replicates -- see
// the README's cost discussion.
func buildNullDistributions(ctx context.Context, cfg Config, realPool *seqpool.Pool) (map[string]map[string][]float64, error) {
	// nullMeanCross[origin][scoredAgainst] = []float64, one entry per
	// replicate. Pre-allocate outer maps; inner slices grow as replicates
	// complete. Replicates run sequentially (each one's model-scoring
	// step is itself parallelised across models), since shuffle
	// replicates are cheap in Go-side compute and the expensive part
	// (cmscore calls) is already parallelised within each replicate.
	nullMeanCross := make(map[string]map[string][]float64)
	origins := realPool.Origins()
	for _, o := range origins {
		nullMeanCross[o] = make(map[string][]float64)
		for _, m := range cfg.Models {
			nullMeanCross[o][m.Name] = make([]float64, 0, cfg.NShuffles)
		}
	}

	nullTick := newTicker(cfg.OnProgress, "scoring_null", len(cfg.Models)*cfg.NShuffles)

	for rep := 0; rep < cfg.NShuffles; rep++ {
		repSeed := cfg.Seed + 7_000_003 + int64(rep)*1_000_033
		shuffledPool, err := shuffleWholePool(realPool, repSeed)
		if err != nil {
			return nil, fmt.Errorf("shuffle replicate %d: %w", rep, err)
		}

		repMatrix, err := scorePoolAgainstAllModels(ctx, cfg, shuffledPool, nullTick)
		if err != nil {
			return nil, fmt.Errorf("shuffle replicate %d: %w", rep, err)
		}

		for _, o := range origins {
			idsO := idsOf(shuffledPool.ByOrigin(o))
			for _, m := range cfg.Models {
				mean, err := repMatrix.MeanFor(idsO, m.Name)
				if err != nil {
					return nil, fmt.Errorf("shuffle replicate %d: %w", rep, err)
				}
				nullMeanCross[o][m.Name] = append(nullMeanCross[o][m.Name], mean)
			}
		}
	}
	return nullMeanCross, nil
}

// shuffleWholePool dinucleotide-shuffles every sequence in pool
// independently, preserving IDs and Origins (so downstream lookups by ID
// still work against the shuffled pool exactly as they do against the
// real one) but randomising residues.
func shuffleWholePool(pool *seqpool.Pool, seed int64) (*seqpool.Pool, error) {
	rng := rand.New(rand.NewSource(seed))
	out := &seqpool.Pool{}
	for _, s := range pool.Sequences {
		shuffled, err := shuffle.DinucleotideShuffle(s.Residues, rng)
		if err != nil {
			return nil, fmt.Errorf("shuffling sequence %q: %w", s.ID, err)
		}
		out.Add(seqpool.Sequence{ID: s.ID, Origin: s.Origin, Residues: shuffled})
	}
	return out, nil
}

// computeAllPairStats slices the real matrix and null distributions into
// one stats.PairStat per reported pair.
func computeAllPairStats(cfg Config, realPool *seqpool.Pool, realMatrix *matrix.ScoreMatrix, nullMeanCross map[string]map[string][]float64) ([]stats.PairStat, error) {
	var results []stats.PairStat
	for i := 0; i < len(cfg.Models); i++ {
		for j := i + 1; j < len(cfg.Models); j++ {
			a, b := cfg.Models[i].Name, cfg.Models[j].Name
			if len(cfg.Query) > 0 && !cfg.Query[a] && !cfg.Query[b] {
				continue // neither side of this pair is a queried model
			}
			nullAonB := nullMeanCross[a][b]
			nullBonA := nullMeanCross[b][a]
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

// Command ksmwcompare is cm-simile's third comparison mode: instead of
// bitscoresumdifference's bootstrap-resampled-mean null (reservoir.go)
// or minmaxcmcompare's evolutionary worst-case search, this draws ONE
// real pool and ONE dinucleotide-shuffled null pool per model (same
// size, single instantiation) and compares their raw score distributions
// directly against every other model via two-sample Kolmogorov-Smirnov
// and Mann-Whitney U tests (src/nptest).
//
// Why a single null pool is a reasonable simplification, not just a
// shortcut: for a two-sample test, power is governed by the effective
// sample size n1*n2/(n1+n2), which is bottlenecked by whichever sample
// is smaller. Once a null sample is several times larger than the real
// pool, growing it further buys negligible extra power -- so
// bitscoresumdifference's adaptive reservoir was never buying detection
// power beyond what -npool itself allows; it was buying a more STABLE
// estimate of where the null sits, averaging out one particular shuffle
// instantiation's randomness. A single null pool trades a bit of
// run-to-run reproducibility for a much simpler system with no growth
// loop, no convergence diagnostic, no adaptive sizing -- worth it if you
// don't need bit-for-bit-stable p-values across reruns, or if a few
// independent null pools (not the full reservoir machinery) are enough
// insurance against an unlucky one-off shuffle.
//
// KS and Mann-Whitney test different things and are both reported,
// deliberately not collapsed into one number:
//   - Mann-Whitney: is the real pool typically higher-scoring than the
//     null (a location shift)? Matches the same "does real structure
//     beat noise" hypothesis bitscoresumdifference/reservoir.go tests,
//     just via direct raw-score comparison rather than resampled means.
//   - Kolmogorov-Smirnov: is there ANY distributional difference,
//     including shape? Catches a real pool that's a MIXTURE -- some
//     sequences retaining strong family-specific signal, others looking
//     noise-like (e.g. partial/degraded homology) -- which can produce a
//     bimodal real distribution that a pure location-shift test
//     underpowers. Report D (bounded, sample-size-independent) alongside
//     the p-value, since KS p-values get very sensitive at large sample
//     sizes and can flag distributional wobbles too small to matter
//     biologically.
//
// -combine mode: each round compares every current pair (from the
// shared pool/matrix, see below), flags any pair whose
// max(ks_d_a_into_b, ks_d_b_into_a) is at or above -combine-threshold-ks
// as a merge candidate (the direction with the larger D names the
// parent -- i.e. the model whose real pool the child looks most like),
// greedily picks a disjoint set of candidates (no model merged twice in
// one round) via chooseDisjointKSCandidates, and merges each pair with
// combine.CombineParentWithTopChild. This repeats until a round
// produces no qualifying pairs (models have separated enough) or
// -combine-max-rounds is hit.
//
// The pool and score matrix persist across the whole combine loop
// rather than being rebuilt each round. When a round produces new
// (merged) models, growPoolAndScore only samples sequences for those
// new models and only runs the two score passes needed to fill in
// every entry a comparison could touch: each new model's CM against
// the whole pool (covers new-into-old and new-into-new), and each
// surviving old model's CM against just the newly added sequences
// (covers old-into-new). Old models are never rescored against
// sequences already in the matrix. That keeps a round's cost
// proportional to how many models changed that round, not to the size
// of the whole population.
//
// Build:
//
//	go build -o bin/ksmwcompare src/ksmwcompare.go
//
// Usage:
//
//	ksmwcompare -cm motif17=motif17.cm -cm RF00050=RF00050.cm \
//	            -npool 200 -out results/ksmw.tsv
//
// Usage with combine mode (iteratively merge highly similar models and rerun comparisons until no qualifying pair remains):
//	ksmwcompare -cmfile population.cm \
//	            -npool 200 -combine -combine-threshold-ks 0.8 -combine-top-proportion 0.25 \
//	            -out results/ksmw-combine.tsv
package main

import (
	"context"
	"flag"
	"fmt"
	"math/rand"
	"os"
	"path/filepath"
	"sort"
	"strings"
	"sync"

	"cmsimile/src/combine"
	"cmsimile/src/infernal"
	"cmsimile/src/matrix"
	"cmsimile/src/nptest"
	"cmsimile/src/seqpool"
	"cmsimile/src/shuffle"
)

// cmFlag / cmFileFlag: same small helpers duplicated in each standalone
// program in this directory -- see bitscoresumdifference.go /
// minmaxcmcompare.go for the rationale (compiled as independent files,
// not worth a shared package for this little boilerplate).
type modelSpec struct{ Name, Path string }

type cmFlag []modelSpec

func (c *cmFlag) String() string { return fmt.Sprintf("%v", *c) }
func (c *cmFlag) Set(value string) error {
	name, path, ok := strings.Cut(value, "=")
	if !ok || name == "" || path == "" {
		return fmt.Errorf("expected -cm NAME=PATH, got %q", value)
	}
	*c = append(*c, modelSpec{Name: name, Path: path})
	return nil
}

type cmFileFlag []string

func (c *cmFileFlag) String() string { return strings.Join(*c, ",") }
func (c *cmFileFlag) Set(value string) error {
	*c = append(*c, value)
	return nil
}

func main() {
	if err := run(); err != nil {
		fmt.Fprintln(os.Stderr, "ksmwcompare:", err)
		os.Exit(1)
	}
}

func run() error {
	var models cmFlag
	var cmFiles cmFileFlag
	flag.Var(&models, "cm", "single-model covariance model as NAME=PATH; repeat for each model")
	flag.Var(&cmFiles, "cm-file", "path to a (possibly multi-model) CM file; every named model inside is expanded automatically")
	query := flag.String("query", "", "comma-separated model names to restrict reported pairs to (default: report every pair)")

	nPool := flag.Int("npool", 200, "size of both the real pool and the null (shuffled) pool per model")
	nNullPools := flag.Int("n-null-pools", 1, "independent shuffled null pool instantiations per model, pooled together before testing -- cheap insurance against an unlucky single shuffle without the full reservoir machinery (see package doc comment)")
	seed := flag.Int64("seed", 42, "base random seed (deterministic reproducibility across runs)")
	threads := flag.Int("threads", 4, "max concurrent cmemit/cmsearch subprocess calls dispatched at once")
	outPath := flag.String("out", "", "output TSV path (default: stdout)")

	cmemitBin := flag.String("cmemit-bin", "", "path to cmemit binary (default: resolve \"cmemit\" from PATH)")
	cmsearchBin := flag.String("cmsearch-bin", "", "path to cmsearch binary (default: resolve \"cmsearch\" from PATH)")
	cmstatBin := flag.String("cmstat-bin", "", "path to cmstat binary, used for -cm-file expansion")
	cmfetchBin := flag.String("cmfetch-bin", "", "path to cmfetch binary, used for -cm-file expansion")
	cmsearchCPU := flag.Int("cmsearch-cpu", 1, "threads passed to cmsearch's own --cpu")
	minBitScore := flag.Float64("min-bitscore", -1000.0, "cmsearch -T inclusion threshold; keep low so weak/null sequences are never excluded from scoring")
	noFilters := flag.Bool("disable-filters", true, "pass --max to cmsearch, disabling acceleration filters (required so weak/shuffled sequences always get a real score)")
	keepTemp := flag.Bool("keep-temp", false, "keep temporary FASTA/tblout files instead of deleting them")
	combineFlag := flag.Bool("combine", false, "iteratively combine highly similar model pairs and rerun comparisons until no qualifying pair remains")
	combineThresholdKS := flag.Float64("combine-threshold-ks", 0.8, "combine trigger threshold on max(ks_d_a_into_b, ks_d_b_into_a); pairs at or above this are merge candidates")
	combineTopProportion := flag.Float64("combine-top-proportion", 0.25, "fraction (0,1] of child real-pool sequences (ranked by score vs chosen parent) used when rebuilding the parent")
	combineMaxRounds := flag.Int("combine-max-rounds", 100, "safety cap on combine rounds to prevent accidental non-termination")
	combineOutDir := flag.String("combine-outdir", "", "Output of any combined models is written to this directory (default: working directory)")

	flag.Parse()

	workDir, err := os.MkdirTemp("", "ksmwcompare-*")
	if err != nil {
		return fmt.Errorf("creating working directory: %w", err)
	}
	if !*keepTemp {
		defer os.RemoveAll(workDir)
	} else {
		fmt.Fprintln(os.Stderr, "ksmwcompare: working directory preserved at:", workDir)
	}

	ctx := context.Background()

	expander := infernal.Expander{CmstatBin: *cmstatBin, CmfetchBin: *cmfetchBin}
	for _, cmFile := range cmFiles {
		expanded, err := expander.ExpandCMFile(ctx, cmFile, workDir)
		if err != nil {
			return fmt.Errorf("expanding -cm-file %s: %w", cmFile, err)
		}
		for _, m := range expanded {
			models = append(models, modelSpec{Name: m.Name, Path: m.Path})
		}
	}

	if len(models) < 2 {
		flag.Usage()
		return fmt.Errorf("need at least 2 models total (via -cm and/or -cm-file), got %d", len(models))
	}

	var queryFilter map[string]bool
	if *query != "" {
		if *combineFlag {
			return fmt.Errorf("-combine currently does not support -query; combine mode evaluates and mutates the full model set each round")
		}
		queryFilter = make(map[string]bool)
		for _, name := range strings.Split(*query, ",") {
			name = strings.TrimSpace(name)
			if name != "" {
				queryFilter[name] = true
			}
		}
	}

	emitter := infernal.Emitter{BinPath: *cmemitBin}
	scorer := infernal.Scorer{
		BinPath:        *cmsearchBin,
		CPU:            *cmsearchCPU,
		MinBitScore:    *minBitScore,
		DisableFilters: *noFilters,
		KeepTempFiles:  *keepTemp,
	}

	var rows []resultRow
	if *combineFlag {
		rows, err = runKSMWCombineLoop(ctx, models, *nPool, *nNullPools, *seed, *threads, emitter, scorer, *combineThresholdKS, *combineTopProportion, *combineMaxRounds, combineOutDir)
		if err != nil {
			return err
		}
	} else {
		rows, _, _, err = runKSMWPass(ctx, models, queryFilter, *nPool, *nNullPools, *seed, *threads, emitter, scorer)
		if err != nil {
			return err
		}
	}

	out := os.Stdout
	if *outPath != "" {
		f, err := os.Create(*outPath)
		if err != nil {
			return fmt.Errorf("creating output file: %w", err)
		}
		defer f.Close()
		out = f
	}
	return writeResultsTSV(out, rows)
}

type ksMergeCandidate struct {
	Row       resultRow
	Parent    string
	Child     string
	KSMaxBoth float64
}

func runKSMWCombineLoop(
	ctx context.Context,
	models []modelSpec,
	nPool int,
	nNullPools int,
	seed int64,
	threads int,
	emitter infernal.Emitter,
	scorer infernal.Scorer,
	thresholdKS float64,
	topProportion float64,
	maxRounds int,
	combineOutDir *string,
) ([]resultRow, error) {
	if topProportion <= 0 || topProportion > 1 {
		return nil, fmt.Errorf("-combine-top-proportion must be in (0,1], got %.4f", topProportion)
	}
	if maxRounds <= 0 {
		return nil, fmt.Errorf("-combine-max-rounds must be > 0, got %d", maxRounds)
	}

	workDir := *combineOutDir
	if workDir == "" {
		workDir = "."
	}

	current := append([]modelSpec(nil), models...)
	pool := &seqpool.Pool{}
	scoreMatrix := matrix.New()

	// Seed the pool/matrix for the starting model set: with no old
	// models yet, this is just "sample and fully score every model",
	// same as a single runKSMWPass would do.
	if err := growPoolAndScore(ctx, nil, current, pool, scoreMatrix, nPool, nNullPools, seed, threads, emitter, scorer); err != nil {
		return nil, fmt.Errorf("initial comparison pass: %w", err)
	}

	var lastRows []resultRow

	for round := 1; round <= maxRounds; round++ {
		if len(current) < 2 {
			break
		}

		rows, err := compareAllPairs(pool, scoreMatrix, current, nil)
		if err != nil {
			return nil, fmt.Errorf("combine round %d: comparison: %w", round, err)
		}
		lastRows = rows

		candidates := make([]ksMergeCandidate, 0)
		for _, r := range rows {
			ksMax := r.KSDAintoB
			parent, child := r.ModelB, r.ModelA
			if r.KSDBintoA > ksMax {
				ksMax = r.KSDBintoA
				parent, child = r.ModelA, r.ModelB
			}
			if ksMax >= thresholdKS {
				candidates = append(candidates, ksMergeCandidate{
					Row:       r,
					Parent:    parent,
					Child:     child,
					KSMaxBoth: ksMax,
				})
			}
		}

		if len(candidates) == 0 {
			// No pairs qualify for merging this round -- the population
			// has separated enough. Stop iterating.
			break
		}

		// Merge the most-similar pairs first, skipping any candidate
		// that reuses a model already claimed by a merge this round.
		sort.Slice(candidates, func(i, j int) bool {
			return candidates[i].KSMaxBoth > candidates[j].KSMaxBoth
		})
		// selected := chooseDisjointKSCandidates(candidates)
		selected := chooseTopKSCandidate(candidates)

		newModels, removedModels, err := applyKSMergeRound(
			ctx,
			current,
			selected,
			pool,
			scoreMatrix,
			workDir,
			topProportion,
			round,
			threads,
		)
		if err != nil {
			return nil, fmt.Errorf("combine round %d: %w", round, err)
		}

		removedSet := make(map[string]bool, len(removedModels))
		for _, m := range removedModels {
			removedSet[m.Name] = true
		}

		survivingOld := make([]modelSpec, 0, len(current)-len(removedModels))
		for _, m := range current {
			if !removedSet[m.Name] {
				survivingOld = append(survivingOld, m)
			}
		}

		// Distinct seed per round so repeated merged-model names across
		// rounds (unlikely, but names aren't guaranteed unique over
		// time) don't resample identically.
		roundSeed := seed + int64(round)*9_700_001
		if err := growPoolAndScore(ctx, survivingOld, newModels, pool, scoreMatrix, nPool, nNullPools, roundSeed, threads, emitter, scorer); err != nil {
			return nil, fmt.Errorf("combine round %d: scoring merged models: %w", round, err)
		}

		current = append(survivingOld, newModels...)
	}

	if len(current) < 2 {
		// Nothing left to compare -- report the last round that had a
		// full set of comparisons.
		return lastRows, nil
	}

	// The last round's `lastRows` were computed before that round's
	// merge, so they're stale relative to the final `current`. This is
	// pure map lookups against the already-populated matrix, no
	// resampling or rescoring.
	finalRows, err := compareAllPairs(pool, scoreMatrix, current, nil)
	if err != nil {
		return nil, fmt.Errorf("final comparison pass: %w", err)
	}
	return finalRows, nil
}

// growPoolAndScore samples real+null sequence pools for newModels, adds
// them to pool, and fills in every score-matrix entry a comparison
// between any two of oldModels∪newModels could need:
//
//   - each new model's CM scored against the WHOLE pool (old+new
//     sequences) -- covers both directions of every new-vs-new pair,
//     and the new-into-old direction of every new-vs-old pair, since
//     every new model does its own full-pool pass.
//   - each old model's CM scored against ONLY the newly added
//     sequences -- covers old-into-new, the one direction a new
//     model's own full-pool pass can't supply.
//
// oldModels' CMs are never rescored against sequences already in the
// pool -- those score-matrix entries exist already from an earlier
// call. That's what keeps a round's cost proportional to how many
// models changed, not to the size of the whole population.
func growPoolAndScore(
	ctx context.Context,
	oldModels []modelSpec,
	newModels []modelSpec,
	pool *seqpool.Pool,
	scoreMatrix *matrix.ScoreMatrix,
	nPool int,
	nNullPools int,
	seed int64,
	threads int,
	emitter infernal.Emitter,
	scorer infernal.Scorer,
) error {
	if len(newModels) == 0 {
		return nil
	}

	newSeqs := &seqpool.Pool{}
	var newSeqsMu sync.Mutex

	if err := parallelFor(len(newModels), threads, func(i int) error {
		m := newModels[i]
		subSeed := seed + int64(i)*1_000_003
		real, err := emitter.Sample(ctx, m.Path, m.Name, nPool, subSeed)
		if err != nil {
			return fmt.Errorf("sampling real pool for %q: %w", m.Name, err)
		}

		nullOrigin := m.Name + "__null"
		var null []seqpool.Sequence
		rng := rand.New(rand.NewSource(seed + 8_000_009 + int64(i)*97))
		for rep := 0; rep < nNullPools; rep++ {
			for k, s := range real {
				shuffled, err := shuffle.DinucleotideShuffle(s.Residues, rng)
				if err != nil {
					return fmt.Errorf("shuffling null sequence for %q: %w", m.Name, err)
				}
				id := fmt.Sprintf("%s_null_%02d_%04d", m.Name, rep, k)
				null = append(null, seqpool.NewSequence(id, nullOrigin, shuffled))
			}
		}

		poolMu.Lock()
		pool.Add(real...)
		pool.Add(null...)
		poolMu.Unlock()

		newSeqsMu.Lock()
		newSeqs.Add(real...)
		newSeqs.Add(null...)
		newSeqsMu.Unlock()
		return nil
	}); err != nil {
		return err
	}

	// Each new model scores the FULL pool (old sequences + everything
	// just added). This alone supplies both directions of every
	// new-vs-new pair, and the new-into-old direction of every
	// new-vs-old pair.
	if err := parallelFor(len(newModels), threads, func(i int) error {
		m := newModels[i]
		scores, err := scorer.Score(ctx, m.Path, pool.Sequences)
		if err != nil {
			return fmt.Errorf("scoring pool against %q: %w", m.Name, err)
		}
		matrixMu.Lock()
		scoreMatrix.SetAll(m.Name, scores)
		matrixMu.Unlock()
		return nil
	}); err != nil {
		return err
	}

	// Each already-scored old model only needs to score the newly
	// added sequences -- its scores against every pre-existing
	// sequence are already in scoreMatrix from an earlier call.
	if len(newSeqs.Sequences) == 0 || len(oldModels) == 0 {
		return nil
	}
	if err := parallelFor(len(oldModels), threads, func(i int) error {
		m := oldModels[i]
		scores, err := scorer.Score(ctx, m.Path, newSeqs.Sequences)
		if err != nil {
			return fmt.Errorf("scoring new sequences against %q: %w", m.Name, err)
		}
		matrixMu.Lock()
		scoreMatrix.SetAll(m.Name, scores)
		matrixMu.Unlock()
		return nil
	}); err != nil {
		return err
	}

	return nil
}

// compareAllPairs runs comparePair over every pair in models (subject to
// queryFilter), using scores already present in scoreMatrix. Pure
// lookups -- no sampling or subprocess calls.
func compareAllPairs(pool *seqpool.Pool, scoreMatrix *matrix.ScoreMatrix, models []modelSpec, queryFilter map[string]bool) ([]resultRow, error) {
	rows := make([]resultRow, 0)
	for i := 0; i < len(models); i++ {
		for j := i + 1; j < len(models); j++ {
			a, b := models[i], models[j]
			if len(queryFilter) > 0 && !(queryFilter[a.Name] || queryFilter[b.Name]) {
				continue
			}
			row, err := comparePair(pool, scoreMatrix, a.Name, b.Name)
			if err != nil {
				return nil, fmt.Errorf("pair (%s, %s): %w", a.Name, b.Name, err)
			}
			rows = append(rows, row)
		}
	}
	return rows, nil
}

func applyKSMergeRound(
	ctx context.Context,
	models []modelSpec,
	selected []ksMergeCandidate,
	pool *seqpool.Pool,
	scoreMatrix *matrix.ScoreMatrix,
	workDir string,
	topProportion float64,
	round int,
	threads int,
) ([]modelSpec, []modelSpec, error) {

	byName := make(map[string]modelSpec, len(models))
	for _, m := range models {
		byName[m.Name] = m
	}

	// selected is already disjoint (see chooseDisjointKSCandidates), so
	// the merges in it are independent and safe to run concurrently.
	// Calibration parallelises internally (cmcalibrate's own --cpu), so
	// split this round's -threads budget evenly across however many
	// merges are actually running at once, rather than starving every
	// merge down to 1 cpu (sequential) or oversubscribing the machine
	// (each merge getting the full -threads count).
	concurrency := len(selected)
	if concurrency > threads {
		concurrency = threads
	}
	if concurrency < 1 {
		concurrency = 1
	}
	cpuPerMerge := threads / concurrency
	if cpuPerMerge < 1 {
		cpuPerMerge = 1
	}

	// each merge removes BOTH parent and child; index-assigned by
	// worker below so no mutex is needed across goroutines.
	newModels := make([]modelSpec, len(selected))
	removedModels := make([]modelSpec, 2*len(selected))

	err := parallelFor(len(selected), concurrency, func(i int) error {
		c := selected[i]

		parent, ok := byName[c.Parent]
		if !ok {
			return fmt.Errorf("combine: missing parent model %q", c.Parent)
		}

		child, ok := byName[c.Child]
		if !ok {
			return fmt.Errorf("combine: missing child model %q", c.Child)
		}

		newName := combine.CombinedModelName(parent.Name, child.Name)

		outPath := filepath.Join(
			workDir,
			fmt.Sprintf(
				"combine_ks_round_%03d_%02d_%s.cm",
				round,
				i,
				combine.SafeFileStem(newName),
			),
		)

		if _, err := combine.CombineParentWithTopChild(
			ctx,
			parent.Path,
			outPath,
			parent.Name,
			child.Name,
			pool,
			scoreMatrix,
			topProportion,
			cpuPerMerge,
		); err != nil {
			return fmt.Errorf("combine %s <- %s: %w", parent.Name, child.Name, err)
		}

		newModels[i] = modelSpec{Name: newName, Path: outPath}
		// IMPORTANT: remove BOTH original models
		removedModels[2*i] = parent
		removedModels[2*i+1] = child
		return nil
	})
	if err != nil {
		return nil, nil, err
	}

	return newModels, removedModels, nil
}

func chooseDisjointKSCandidates(candidates []ksMergeCandidate) []ksMergeCandidate {
	used := make(map[string]bool)
	out := make([]ksMergeCandidate, 0, len(candidates))
	for _, c := range candidates {
		if used[c.Parent] || used[c.Child] {
			continue
		}
		used[c.Parent] = true
		used[c.Child] = true
		out = append(out, c)
	}
	return out
}

func chooseTopKSCandidate(candidates []ksMergeCandidate) []ksMergeCandidate {
	if len(candidates) == 0 {
		return nil
	}

	// Pick the single candidate with the highest KSMaxBoth value.
	top := candidates[0]
	for _, c := range candidates[1:] {
		if c.KSMaxBoth > top.KSMaxBoth {
			top = c
		}
	}
	return []ksMergeCandidate{top}
}

// runKSMWPass runs a single one-shot comparison: sample and fully score
// every model (there's no "old" set on a first pass), then compare
// every requested pair. Used by non-combine mode; combine mode uses
// growPoolAndScore/compareAllPairs directly so it can grow the pool and
// matrix incrementally across rounds instead of rebuilding them fresh
// each time.
func runKSMWPass(
	ctx context.Context,
	models []modelSpec,
	queryFilter map[string]bool,
	nPool int,
	nNullPools int,
	seed int64,
	threads int,
	emitter infernal.Emitter,
	scorer infernal.Scorer,
) ([]resultRow, *seqpool.Pool, *matrix.ScoreMatrix, error) {
	pool := &seqpool.Pool{}
	scoreMatrix := matrix.New()

	if err := growPoolAndScore(ctx, nil, models, pool, scoreMatrix, nPool, nNullPools, seed, threads, emitter, scorer); err != nil {
		return nil, nil, nil, err
	}

	rows, err := compareAllPairs(pool, scoreMatrix, models, queryFilter)
	if err != nil {
		return nil, nil, nil, err
	}

	return rows, pool, scoreMatrix, nil
}

var poolMu sync.Mutex
var matrixMu sync.Mutex

type resultRow struct {
	ModelA, ModelB string

	MedianAintoB float64
	MedianBintoA float64

	MedianNullintoA float64
	MedianNullintoB float64

	KSDAintoB    float64
	KSPValAintoB float64
	MWAUCAintoB  float64
	MWPValAintoB float64

	KSDBintoA    float64
	KSPValBintoA float64
	MWAUCBintoA  float64
	MWPValBintoA float64

	N1, N2 int
}

// comparePair compares A's real-vs-null scores against B, and B's
// real-vs-null scores against A, via both tests. Each direction asks:
// "when scored against the OTHER model, does this model's real pool
// look different from its own null pool?" -- i.e. does A's real
// structure earn more recognition from B than composition-matched noise
// would, tested directly on raw scores rather than resampled means.
func comparePair(pool *seqpool.Pool, m *matrix.ScoreMatrix, a, b string) (resultRow, error) {
	realA := scoresFor(pool, m, a, b)
	nullA := scoresFor(pool, m, a+"__null", b)
	realB := scoresFor(pool, m, b, a)
	nullB := scoresFor(pool, m, b+"__null", a)

	if len(realA) == 0 || len(nullA) == 0 || len(realB) == 0 || len(nullB) == 0 {
		return resultRow{}, fmt.Errorf("missing scores for pair (%s, %s) -- expected non-empty real and null pools for both models", a, b)
	}

	ksAB := nptest.TwoSampleKS(realA, nullA)
	mwAB := nptest.MannWhitneyU(realA, nullA)
	ksBA := nptest.TwoSampleKS(realB, nullB)
	mwBA := nptest.MannWhitneyU(realB, nullB)

	return resultRow{
		ModelA: a, ModelB: b,
		MedianAintoB: median(realA), MedianBintoA: median(realB),
		MedianNullintoA: median(nullA), MedianNullintoB: median(nullB),
		KSDAintoB: ksAB.D, KSPValAintoB: ksAB.PValue, MWAUCAintoB: mwAB.AUC, MWPValAintoB: mwAB.PValue,
		KSDBintoA: ksBA.D, KSPValBintoA: ksBA.PValue, MWAUCBintoA: mwBA.AUC, MWPValBintoA: mwBA.PValue,
		N1: len(realA), N2: len(nullA),
	}, nil
}

func scoresFor(pool *seqpool.Pool, m *matrix.ScoreMatrix, origin, target string) []float64 {
	seqs := pool.ByOrigin(origin)
	out := make([]float64, 0, len(seqs))
	for _, s := range seqs {
		if sc, ok := m.Get(s.ID, target); ok {
			out = append(out, sc)
		}
	}
	return out
}

func median(xs []float64) float64 {
	if len(xs) == 0 {
		return 0.0
	}
	sorted := make([]float64, len(xs))
	copy(sorted, xs)
	sort.Float64s(sorted)
	mid := len(sorted) / 2
	if len(sorted)%2 == 0 {
		return (sorted[mid-1] + sorted[mid]) / 2.0
	}
	return sorted[mid]
}

func writeResultsTSV(w *os.File, rows []resultRow) error {
	header := "model_a\tmodel_b\tn_real\tn_null\t" +
		"median_a_into_b\tmedian_b_into_a\tmedian_null_into_a\tmedian_null_into_b\t" +
		"ks_d_a_into_b\tks_pvalue_a_into_b\tmw_auc_a_into_b\tmw_pvalue_a_into_b\t" +
		"ks_d_b_into_a\tks_pvalue_b_into_a\tmw_auc_b_into_a\tmw_pvalue_b_into_a"
	if _, err := fmt.Fprintln(w, header); err != nil {
		return err
	}
	for _, r := range rows {
		if _, err := fmt.Fprintf(w, "%s\t%s\t%d\t%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6g\t%.6f\t%.6g\t%.6f\t%.6g\t%.6f\t%.6g\n",
			r.ModelA, r.ModelB, r.N1, r.N2,
			r.MedianAintoB, r.MedianBintoA, r.MedianNullintoA, r.MedianNullintoB,
			r.KSDAintoB, r.KSPValAintoB, r.MWAUCAintoB, r.MWPValAintoB,
			r.KSDBintoA, r.KSPValBintoA, r.MWAUCBintoA, r.MWPValBintoA,
		); err != nil {
			return err
		}
	}
	return nil
}

// parallelFor: same bounded-concurrency helper duplicated in every
// standalone program in this directory -- see minmaxcmcompare.go.
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
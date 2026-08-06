// Command minmaxcmcompare estimates the empirical MaxiMin similarity
// between Infernal covariance models: for each pair, it searches for the
// single sequence that scores as well as possible against BOTH models
// simultaneously (max over candidates of min(scoreA, scoreB)), the same
// quantity CMCompare (Höner zu Siederdissen & Hofacker 2010) computes
// exactly via dynamic programming over the product of the two CMs' state
// graphs. This tool instead *approximates* it: seed a population of real
// cmemit-sampled sequences, then refine with an evolutionary search
// (mutation + crossover, selection by fitness = min(scoreA, scoreB) via
// cmsearch) toward the worst-case-best sequence.
//
// This is the worst-case counterpart to bitscoresumdifference.go's
// distributional/enrichment statistic in this same directory: maximin
// answers "does even one dangerous sequence exist that would confuse
// both models", not "on average, how much do these models overlap". See
// the top-level README's Algorithm Approach section for when each is the
// right tool. Because a maximin search needs live mutate-then-rescore
// feedback per candidate, it does not batch across many pairs the way
// the bootstrap approach does -- expect this tool to cost meaningfully
// more per pair, and to be used for confirming a small number of
// promising hits (e.g. ones bitscoresumdifference already flagged as
// significant) rather than screening thousands of pairs from scratch.
//
// Like bitscoresumdifference.go, this file is its own standalone
// `package main` living in the same directory as other standalone
// programs -- build it by naming the file explicitly:
//
//	go build -o bin/minmaxcmcompare src/minmaxcmcompare.go
//
// Usage:
//
//	minmaxcmcompare -cm motif17=motif17.cm -cm RF00050=RF00050.cm \
//	                 -generations 20 -pop-size 40 -out results/maximin.tsv
package main

import (
	"context"
	"flag"
	"fmt"
	"math"
	"math/rand"
	"os"
	"sort"
	"strings"

	"cmsimile/src/infernal"
	"cmsimile/src/seqpool"
)

// cmFlag / cmFileFlag: identical small helpers to bitscoresumdifference.go.
// Duplicated rather than shared, since these two programs are compiled as
// independent standalone files (see the package doc comment on both) and
// this boilerplate is small enough that a shared package would add more
// indirection than it saves.
type cmFlag []modelSpec

type modelSpec struct {
	Name string
	Path string
}

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
	if err := runMinmaxcmcompare(); err != nil {
		fmt.Fprintln(os.Stderr, "minmaxcmcompare:", err)
		os.Exit(1)
	}
}

func runMinmaxcmcompare() error {
	var models cmFlag
	var cmFiles cmFileFlag
	flag.Var(&models, "cm", "single-model covariance model as NAME=PATH; repeat for each model")
	flag.Var(&cmFiles, "cm-file", "path to a (possibly multi-model) CM file; every named model inside is expanded automatically")
	query := flag.String("query", "", "comma-separated model names to restrict reported pairs to (default: report every pair)")

	popSize := flag.Int("pop-size", 40, "evolutionary population size per pair")
	generations := flag.Int("generations", 15, "number of mutate/select generations per pair")
	survivalFrac := flag.Float64("survival-frac", 0.25, "fraction of the population retained as survivors/parents each generation")
	nSeed := flag.Int("nseed", 15, "cmemit-sampled sequences drawn per model to seed the initial population, in addition to each model's consensus sequence")
	mutationsPerChild := flag.Int("mutations-per-child", 2, "point mutations applied when producing a child from a parent")
	crossoverProb := flag.Float64("crossover-prob", 0.3, "probability a child is produced by crossing two survivors rather than mutating one")

	seed := flag.Int64("seed", 42, "base random seed (deterministic reproducibility across runs)")
	outPath := flag.String("out", "", "output TSV path (default: stdout)")

	cmemitBin := flag.String("cmemit-bin", "", "path to cmemit binary (default: resolve \"cmemit\" from PATH)")
	cmsearchBin := flag.String("cmsearch-bin", "", "path to cmsearch binary (default: resolve \"cmsearch\" from PATH)")
	cmstatBin := flag.String("cmstat-bin", "", "path to cmstat binary, used for -cm-file expansion")
	cmfetchBin := flag.String("cmfetch-bin", "", "path to cmfetch binary, used for -cm-file expansion")
	cmsearchCPU := flag.Int("cmsearch-cpu", 1, "threads passed to cmsearch's own --cpu")
	minBitScore := flag.Float64("min-bitscore", -1000.0, "cmsearch -T inclusion threshold; keep low so weak candidates are never excluded from scoring")
	noFilters := flag.Bool("disable-filters", true, "pass --max to cmsearch, disabling acceleration filters (slower, but required so weak/mutated candidates always get a real score)")
	keepTemp := flag.Bool("keep-temp", false, "keep temporary FASTA/tblout files instead of deleting them")

	flag.Parse()

	workDir, err := os.MkdirTemp("", "minmaxcmcompare-*")
	if err != nil {
		return fmt.Errorf("creating working directory: %w", err)
	}
	if !*keepTemp {
		defer os.RemoveAll(workDir)
	} else {
		fmt.Fprintln(os.Stderr, "minmaxcmcompare: working directory preserved at:", workDir)
	}

	expander := infernal.Expander{CmstatBin: *cmstatBin, CmfetchBin: *cmfetchBin}
	for _, cmFile := range cmFiles {
		expanded, err := expander.ExpandCMFile(context.Background(), cmFile, workDir)
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

	ctx := context.Background()

	// Seed each model's own pool once, up front, and reuse it as both
	// the evolutionary seed population and the baseline (self-score)
	// reference for every pair that model participates in.
	seedPools := make(map[string][]seqpool.Sequence, len(models))
	selfScores := make(map[string]float64, len(models))
	for i, m := range models {
		subSeed := *seed + int64(i)*1_000_003
		pool, err := seedPool(ctx, emitter, scorer, m, *nSeed, subSeed)
		if err != nil {
			return fmt.Errorf("seeding pool for model %q: %w", m.Name, err)
		}
		seedPools[m.Name] = pool.seqs
		selfScores[m.Name] = pool.meanSelfScore
	}

	var rows []resultRow
	for i := 0; i < len(models); i++ {
		for j := i + 1; j < len(models); j++ {
			a, b := models[i], models[j]
			if len(queryFilter) > 0 && !queryFilter[a.Name] && !queryFilter[b.Name] {
				continue
			}
			fmt.Fprintf(os.Stderr, "minmaxcmcompare: refining %s vs %s ...\n", a.Name, b.Name)

			pairSeed := *seed + int64(i)*7_919 + int64(j)*104_729
			result, err := refineMaximin(ctx, emitter, scorer, a, b,
				seedPools[a.Name], seedPools[b.Name],
				refineConfig{
					PopSize:           *popSize,
					Generations:       *generations,
					SurvivalFrac:      *survivalFrac,
					MutationsPerChild: *mutationsPerChild,
					CrossoverProb:     *crossoverProb,
					Seed:              pairSeed,
				})
			if err != nil {
				return fmt.Errorf("refining pair (%s, %s): %w", a.Name, b.Name, err)
			}

			baseline := (selfScores[a.Name] + selfScores[b.Name]) / 2.0
			similarity := 0.0
			if baseline != 0 {
				similarity = result.BestFitness / baseline
			}

			rows = append(rows, resultRow{
				ModelA: a.Name, ModelB: b.Name,
				EmpiricalMaximinScore: result.BestFitness,
				EmpiricalMaximinSeq:   result.BestSeq,
				Baseline:              baseline,
				SimilarityRatio:       similarity,
				Generations:           *generations,
			})
		}
	}

	sort.Slice(rows, func(i, j int) bool { return rows[i].SimilarityRatio > rows[j].SimilarityRatio })

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

// --- seeding ---------------------------------------------------------

type modelPool struct {
	seqs          []seqpool.Sequence
	meanSelfScore float64
}

// seedPool draws a model's consensus sequence plus nSeed sampled
// sequences, scores them against the model's own CM, and returns both
// the pool and its mean self-score (used later as the similarity-ratio
// baseline, exactly as in the original 1-to-1 MaxiMin design).
func seedPool(ctx context.Context, emitter infernal.Emitter, scorer infernal.Scorer, m modelSpec, nSeed int, seed int64) (modelPool, error) {
	var seqs []seqpool.Sequence

	consensus, err := emitter.Consensus(ctx, m.Path, m.Name)
	if err != nil {
		return modelPool{}, fmt.Errorf("cmemit consensus: %w", err)
	}
	seqs = append(seqs, consensus)

	if nSeed > 0 {
		samples, err := emitter.Sample(ctx, m.Path, m.Name, nSeed, seed)
		if err != nil {
			return modelPool{}, fmt.Errorf("cmemit sample: %w", err)
		}
		seqs = append(seqs, samples...)
	}

	scores, err := scorer.Score(ctx, m.Path, seqs)
	if err != nil {
		return modelPool{}, fmt.Errorf("scoring seed pool against own model: %w", err)
	}
	var sum float64
	for _, s := range seqs {
		sum += scores[s.ID]
	}
	return modelPool{seqs: seqs, meanSelfScore: sum / float64(len(seqs))}, nil
}

// --- evolutionary refinement ------------------------------------------

type refineConfig struct {
	PopSize           int
	Generations       int
	SurvivalFrac      float64
	MutationsPerChild int
	CrossoverProb     float64
	Seed              int64
}

type refineResult struct {
	BestFitness float64
	BestSeq     string
}

// refineMaximin runs the evolutionary search described in the package
// doc comment for one pair of models, returning the best (highest
// min(scoreA, scoreB)) sequence found.
//
// Each generation makes exactly two batched cmsearch calls (one against
// modelA, one against modelB, scoring every not-yet-scored candidate at
// once) regardless of population size, keeping the per-generation
// subprocess cost constant rather than growing with PopSize.
func refineMaximin(ctx context.Context, emitter infernal.Emitter, scorer infernal.Scorer, a, b modelSpec, seedA, seedB []seqpool.Sequence, cfg refineConfig) (refineResult, error) {
	rng := rand.New(rand.NewSource(cfg.Seed))

	type cand struct {
		id             string
		seq            string
		scoreA, scoreB float64
	}

	// Seed the initial population from both models' real sequences --
	// this is the same intuition as the original sampling-only design
	// (typical sequences from either model are a reasonable starting
	// point), refined from there by selection pressure toward sequences
	// good for both.
	var population []cand
	idCounter := 0
	nextID := func() string {
		idCounter++
		return fmt.Sprintf("%s_x_%s_cand_%04d", a.Name, b.Name, idCounter)
	}
	for _, s := range seedA {
		population = append(population, cand{id: nextID(), seq: s.Residues})
	}
	for _, s := range seedB {
		population = append(population, cand{id: nextID(), seq: s.Residues})
	}

	scoreBatch := func(candidates []cand) error {
		if len(candidates) == 0 {
			return nil
		}
		seqs := make([]seqpool.Sequence, len(candidates))
		for i, c := range candidates {
			seqs[i] = seqpool.NewSequence(c.id, "candidate", c.seq)
		}
		scoresA, err := scorer.Score(ctx, a.Path, seqs)
		if err != nil {
			return fmt.Errorf("scoring candidates against %s: %w", a.Name, err)
		}
		scoresB, err := scorer.Score(ctx, b.Path, seqs)
		if err != nil {
			return fmt.Errorf("scoring candidates against %s: %w", b.Name, err)
		}
		for i := range candidates {
			candidates[i].scoreA = scoresA[candidates[i].id]
			candidates[i].scoreB = scoresB[candidates[i].id]
		}
		return nil
	}

	if err := scoreBatch(population); err != nil {
		return refineResult{}, err
	}

	fitness := func(c cand) float64 { return math.Min(c.scoreA, c.scoreB) }

	var best cand
	best.scoreA, best.scoreB = math.Inf(-1), math.Inf(-1)
	updateBest := func() {
		for _, c := range population {
			if fitness(c) > fitness(best) {
				best = c
			}
		}
	}
	updateBest()

	survivorCount := int(math.Ceil(float64(cfg.PopSize) * cfg.SurvivalFrac))
	if survivorCount < 1 {
		survivorCount = 1
	}

	for gen := 0; gen < cfg.Generations; gen++ {
		sort.Slice(population, func(i, j int) bool { return fitness(population[i]) > fitness(population[j]) })
		if len(population) > survivorCount {
			population = population[:survivorCount]
		}

		var children []cand
		for len(children) < cfg.PopSize-len(population) {
			var childSeq string
			if len(population) >= 2 && rng.Float64() < cfg.CrossoverProb {
				p1 := population[rng.Intn(len(population))]
				p2 := population[rng.Intn(len(population))]
				childSeq = crossover(p1.seq, p2.seq, rng)
			} else {
				parent := population[rng.Intn(len(population))]
				childSeq = mutate(parent.seq, rng, cfg.MutationsPerChild)
			}
			children = append(children, cand{id: nextID(), seq: childSeq})
		}

		if err := scoreBatch(children); err != nil {
			return refineResult{}, fmt.Errorf("generation %d: %w", gen, err)
		}
		population = append(population, children...)
		updateBest()
	}

	return refineResult{BestFitness: fitness(best), BestSeq: best.seq}, nil
}

var alphabet = []byte("ACGT")

// mutate returns a copy of seq with n random point substitutions.
func mutate(seq string, rng *rand.Rand, n int) string {
	if len(seq) == 0 {
		return seq
	}
	b := []byte(seq)
	for i := 0; i < n; i++ {
		pos := rng.Intn(len(b))
		b[pos] = alphabet[rng.Intn(len(alphabet))]
	}
	return string(b)
}

// crossover splices two sequences at independently-chosen proportional
// split points (rather than requiring equal length, since candidates
// descended from different models' pools may differ in length -- see the
// package doc comment on why this doesn't need positional correspondence
// between the two parents).
func crossover(a, b string, rng *rand.Rand) string {
	if len(a) == 0 {
		return b
	}
	if len(b) == 0 {
		return a
	}
	frac := rng.Float64()
	splitA := int(frac * float64(len(a)))
	splitB := int(frac * float64(len(b)))
	return a[:splitA] + b[splitB:]
}

// --- output ------------------------------------------------------------

type resultRow struct {
	ModelA, ModelB        string
	EmpiricalMaximinScore float64
	EmpiricalMaximinSeq   string
	Baseline              float64
	SimilarityRatio       float64
	Generations           int
}

func writeResultsTSV(w *os.File, rows []resultRow) error {
	if _, err := fmt.Fprintln(w, "model_a\tmodel_b\tsimilarity_ratio\tempirical_maximin_score\tempirical_maximin_seq\tbaseline\tgenerations"); err != nil {
		return err
	}
	for _, r := range rows {
		if _, err := fmt.Fprintf(w, "%s\t%s\t%.6f\t%.6f\t%s\t%.6f\t%d\n",
			r.ModelA, r.ModelB, r.SimilarityRatio, r.EmpiricalMaximinScore, r.EmpiricalMaximinSeq, r.Baseline, r.Generations); err != nil {
			return err
		}
	}
	return nil
}

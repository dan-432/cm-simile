// Command bitscoresumdifference estimates sequence/structure similarity
// between Infernal covariance models via directional mean bit-score
// difference, tested against a dinucleotide-shuffled null distribution
// for a defensible bootstrap p-value per pair. This is cm-simile's
// primary, fast screening tool -- see minmaxcmcompare.go for the
// alternative worst-case (empirical MaxiMin) comparison in this same
// directory.
//
// Because this file and minmaxcmcompare.go are both `package main` with
// their own `func main`, they are NOT built together as a whole-directory
// package -- build (or vet/test) each standalone program by naming its
// file explicitly:
//
//	go build -o bin/bitscoresumdifference src/bitscoresumdifference.go
//
// (`go build ./src/...` will fail with "main redeclared" precisely
// because both programs live in the same directory; this is expected --
// see the top-level README for the reasoning and the full build
// commands for both tools.)
//
// Usage, single-model CM files:
//
//	bitscoresumdifference -cm motif17=path/to/motif17.cm -cm RF00050=path/to/RF00050.cm \
//	         -npool 200 -nshuffles 200 -out results.tsv
//
// Usage, comparing against a multi-model file directly (e.g. Rfam.cm) --
// every named model inside it is expanded into its own model
// automatically, no pre-splitting required:
//
//	bitscoresumdifference -cm motif17=path/to/motif17.cm -cm-file data/Rfam.cm \
//	         -query motif17 -out results/screen.tsv
//
// -query restricts which *pairs* are written to output, but every
// supplied/expanded model is still pooled and scored -- see the
// pipeline package doc comment for the cost model.
package main

import (
	"context"
	"flag"
	"fmt"
	"os"
	"path/filepath"
	"sort"
	"strings"
	"sync"

	"cmsimile/src/combine"
	"cmsimile/src/infernal"
	"cmsimile/src/matrix"
	"cmsimile/src/output"
	"cmsimile/src/pipeline"
	"cmsimile/src/seqpool"
	"cmsimile/src/stats"
)

// cmFlag accumulates repeated -cm name=path flags into a slice.
type cmFlag []pipeline.Model

func (c *cmFlag) String() string { return fmt.Sprintf("%v", *c) }

func (c *cmFlag) Set(value string) error {
	name, path, ok := strings.Cut(value, "=")
	if !ok || name == "" || path == "" {
		return fmt.Errorf("expected -cm NAME=PATH, got %q", value)
	}
	*c = append(*c, pipeline.Model{Name: name, Path: path})
	return nil
}

// cmFileFlag accumulates repeated -cm-file path flags into a slice.
type cmFileFlag []string

func (c *cmFileFlag) String() string { return strings.Join(*c, ",") }
func (c *cmFileFlag) Set(value string) error {
	*c = append(*c, value)
	return nil
}

func main() {
	if err := runBitscoresumdifference(); err != nil {
		fmt.Fprintln(os.Stderr, "bitscoresumdifference:", err)
		os.Exit(1)
	}
}

func runBitscoresumdifference() error {
	var models cmFlag
	var cmFiles cmFileFlag
	flag.Var(&models, "cm", "single-model covariance model as NAME=PATH; repeat for each model")
	flag.Var(&cmFiles, "cm-file", "path to a (possibly multi-model) CM file, e.g. data/Rfam.cm; every named model inside is expanded automatically. Repeat to combine multiple files")
	query := flag.String("query", "", "comma-separated model names to restrict reported pairs to (default: report every pair)")
	nPool := flag.Int("npool", 200, "real sequences sampled per model via cmemit")
	nShuffles := flag.Int("nshuffles", 200, "bootstrap resamples per pair, drawn in-memory from the null reservoir (cheap -- see -reservoir-* flags; no longer drives cmsearch call count the way it used to)")
	resampleSize := flag.Int("resample-size", 0, "size of each bootstrap draw from the null reservoir; defaults to -npool if 0, keeping null and real means on the same sample size")
	reservoirMin := flag.Int("reservoir-min", 500, "starting size of the null reservoir per model (also the floor -resample-size must not exceed)")
	reservoirMax := flag.Int("reservoir-max", 8000, "cap on reservoir growth; a model hitting this without its convergence diagnostic stabilizing gets a printed warning rather than silently proceeding")
	reservoirTolerance := flag.Float64("reservoir-tolerance", 0.02, "max relative change in the convergence diagnostic between doubling rounds that counts as converged")
	diagnosticResamples := flag.Int("diagnostic-resamples", 200, "bootstrap resamples used only to test reservoir convergence each round (separate from, usually smaller than, -nshuffles)")
	seed := flag.Int64("seed", 42, "base random seed (deterministic reproducibility across runs)")
	threads := flag.Int("threads", 4, "max concurrent cmemit/cmsearch subprocess calls dispatched at once")
	outPath := flag.String("out", "", "output TSV path (default: stdout)")
	dumpMatricesDir := flag.String("dump-matrices", "", "if set, write the full real-pool and null-reservoir (sequence x model) score matrices as long-format TSVs to this directory -- e.g. for plotting the null distribution's actual shape")

	cmemitBin := flag.String("cmemit-bin", "", "path to cmemit binary (default: resolve \"cmemit\" from PATH)")
	cmsearchBin := flag.String("cmsearch-bin", "", "path to cmsearch binary (default: resolve \"cmsearch\" from PATH)")
	cmstatBin := flag.String("cmstat-bin", "", "path to cmstat binary, used for -cm-file expansion (default: resolve \"cmstat\" from PATH)")
	cmfetchBin := flag.String("cmfetch-bin", "", "path to cmfetch binary, used for -cm-file expansion (default: resolve \"cmfetch\" from PATH)")

	cmsearchCPU := flag.Int("cmsearch-cpu", 1, "threads passed to cmsearch's own --cpu (internal to a single cmsearch call; multiplies with -threads)")
	minBitScore := flag.Float64("min-bitscore", -1000.0, "cmsearch -T inclusion threshold; must stay low enough that even weak/shuffled sequences are never excluded from scoring (see Scorer doc comment)")
	noFilters := flag.Bool("disable-filters", true, "pass --max to cmsearch, disabling acceleration filters so every sequence is fully scored by the CM (slower, but required for correct null-distribution statistics; see Scorer doc comment)")

	keepTemp := flag.Bool("keep-temp", false, "keep temporary FASTA/tblout files instead of deleting them (for debugging a parsing mismatch against real Infernal output)")
	quiet := flag.Bool("quiet", false, "suppress the run overview and progress bar (written to stderr; stdout/-out results are unaffected either way)")
	combineFlag := flag.Bool("combine", false, "iteratively combine strongly overlapping model pairs and rerun comparisons until no qualifying pair remains")
	combineThresholdMedian := flag.Float64("combine-threshold-median", 0.0, "combine trigger threshold on min(median_a_on_b, median_b_on_a); pairs at or above this are merge candidates")
	combineTopProportion := flag.Float64("combine-top-proportion", 0.25, "fraction (0,1] of child real-pool sequences (ranked by score vs chosen parent) used when rebuilding the parent")
	combineMaxRounds := flag.Int("combine-max-rounds", 100, "safety cap on combine rounds to prevent accidental non-termination")

	flag.Parse()

	workDir, err := os.MkdirTemp("", "bitscoresumdifference-*")
	if err != nil {
		return fmt.Errorf("creating working directory: %w", err)
	}
	if !*keepTemp {
		defer os.RemoveAll(workDir)
	} else {
		fmt.Fprintln(os.Stderr, "bitscoresumdifference: working directory preserved at:", workDir)
	}

	expander := infernal.Expander{CmstatBin: *cmstatBin, CmfetchBin: *cmfetchBin}
	for _, cmFile := range cmFiles {
		expanded, err := expander.ExpandCMFile(context.Background(), cmFile, workDir)
		if err != nil {
			return fmt.Errorf("expanding -cm-file %s: %w", cmFile, err)
		}
		for _, m := range expanded {
			models = append(models, pipeline.Model{Name: m.Name, Path: m.Path})
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
			if name == "" {
				continue
			}
			queryFilter[name] = true
		}
	}

	cfg := pipeline.Config{
		Models:    models,
		Query:     queryFilter,
		NPool:     *nPool,
		NShuffles: *nShuffles,
		Reservoir: pipeline.ReservoirConfig{
			ResampleSize:        *resampleSize,
			ReservoirMin:        *reservoirMin,
			ReservoirMax:        *reservoirMax,
			ReservoirTolerance:  *reservoirTolerance,
			DiagnosticResamples: *diagnosticResamples,
		},
		Seed:            *seed,
		Threads:         *threads,
		DumpMatricesDir: *dumpMatricesDir,
		Emitter:         infernal.Emitter{BinPath: *cmemitBin},
		Scorer: infernal.Scorer{
			BinPath:        *cmsearchBin,
			CPU:            *cmsearchCPU,
			MinBitScore:    *minBitScore,
			DisableFilters: *noFilters,
			KeepTempFiles:  *keepTemp,
		},
	}

	if !*quiet {
		printOverview(os.Stderr, cfg)
		pb := newProgressBar(os.Stderr)
		cfg.OnProgress = pb.onProgress
	}

	ctx := context.Background()
	var results []stats.PairStat
	var warnings []string
	if *combineFlag {
		results, warnings, err = runBitscoreCombineLoop(ctx, cfg, workDir, *combineThresholdMedian, *combineTopProportion, *combineMaxRounds, *quiet)
		if err != nil {
			return err
		}
	} else {
		results, warnings, err = pipeline.Run(ctx, cfg)
		if err != nil {
			return err
		}
	}
	for _, w := range warnings {
		fmt.Fprintln(os.Stderr, "bitscoresumdifference: warning:", w)
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
	return output.WriteTSV(out, results)
}

type medianMergeCandidate struct {
	Pair           stats.PairStat
	ParentName     string
	ChildName      string
	ThresholdScore float64
}

func runBitscoreCombineLoop(
	ctx context.Context,
	baseCfg pipeline.Config,
	workDir string,
	thresholdMedian float64,
	topProportion float64,
	maxRounds int,
	quiet bool,
) ([]stats.PairStat, []string, error) {
	if topProportion <= 0 || topProportion > 1 {
		return nil, nil, fmt.Errorf("-combine-top-proportion must be in (0,1], got %.4f", topProportion)
	}
	if maxRounds <= 0 {
		return nil, nil, fmt.Errorf("-combine-max-rounds must be > 0, got %d", maxRounds)
	}

	models := append([]pipeline.Model(nil), baseCfg.Models...)
	var allWarnings []string
	var lastResults []stats.PairStat

	for round := 1; round <= maxRounds; round++ {
		if len(models) < 2 {
			break
		}

		cfg := baseCfg
		cfg.Models = models
		cfg.Query = nil

		detailed, err := pipeline.RunDetailed(ctx, cfg)
		if err != nil {
			return nil, allWarnings, fmt.Errorf("combine round %d: %w", round, err)
		}
		allWarnings = append(allWarnings, detailed.Warnings...)
		lastResults = detailed.Results

		candidates := selectMedianMergeCandidates(detailed.Results, thresholdMedian)
		if len(candidates) == 0 {
			if !quiet {
				fmt.Fprintf(os.Stderr, "bitscoresumdifference: combine round %d: no pairs met threshold %.4f; stopping\n", round, thresholdMedian)
			}
			break
		}

		selected := chooseDisjointMedianCandidates(candidates)
		if len(selected) == 0 {
			break
		}

		if !quiet {
			fmt.Fprintf(os.Stderr, "bitscoresumdifference: combine round %d: merging %d pair(s)\n", round, len(selected))
		}

		next, err := applyMedianMergeRound(ctx, models, selected, detailed.RealPool, detailed.RealMatrix, workDir, topProportion, round)
		if err != nil {
			return nil, allWarnings, fmt.Errorf("combine round %d: %w", round, err)
		}
		models = next
	}

	if len(models) < 2 {
		if !quiet {
			fmt.Fprintln(os.Stderr, "bitscoresumdifference: combine mode collapsed model set to <2 models; no final pairwise output")
		}
		return nil, allWarnings, nil
	}

	finalCfg := baseCfg
	finalCfg.Models = models
	results, warnings, err := pipeline.Run(ctx, finalCfg)
	if err != nil {
		return nil, allWarnings, fmt.Errorf("final post-combine run: %w", err)
	}
	allWarnings = append(allWarnings, warnings...)
	if len(results) == 0 {
		results = lastResults
	}
	return results, allWarnings, nil
}

func selectMedianMergeCandidates(results []stats.PairStat, threshold float64) []medianMergeCandidate {
	out := make([]medianMergeCandidate, 0)
	for _, r := range results {
		score := r.MedianScoreAonB
		if r.MedianScoreBonA < score {
			score = r.MedianScoreBonA
		}
		if score < threshold {
			continue
		}

		parent := r.ModelA
		child := r.ModelB
		if r.StatAintoB >= r.StatBintoA {
			parent = r.ModelB
			child = r.ModelA
		}

		out = append(out, medianMergeCandidate{
			Pair:           r,
			ParentName:     parent,
			ChildName:      child,
			ThresholdScore: score,
		})
	}

	sort.Slice(out, func(i, j int) bool {
		if out[i].ThresholdScore != out[j].ThresholdScore {
			return out[i].ThresholdScore > out[j].ThresholdScore
		}
		if out[i].ParentName != out[j].ParentName {
			return out[i].ParentName < out[j].ParentName
		}
		return out[i].ChildName < out[j].ChildName
	})
	return out
}

func chooseDisjointMedianCandidates(candidates []medianMergeCandidate) []medianMergeCandidate {
	used := make(map[string]bool)
	out := make([]medianMergeCandidate, 0, len(candidates))
	for _, c := range candidates {
		if used[c.ParentName] || used[c.ChildName] {
			continue
		}
		used[c.ParentName] = true
		used[c.ChildName] = true
		out = append(out, c)
	}
	return out
}

func applyMedianMergeRound(
	ctx context.Context,
	models []pipeline.Model,
	selected []medianMergeCandidate,
	realPool *seqpool.Pool,
	realMatrix *matrix.ScoreMatrix,
	workDir string,
	topProportion float64,
	round int,
) ([]pipeline.Model, error) {
	byName := make(map[string]pipeline.Model, len(models))
	for _, m := range models {
		byName[m.Name] = m
	}

	replaced := make(map[string]pipeline.Model)
	removed := make(map[string]bool)
	for i, sel := range selected {
		parent, ok := byName[sel.ParentName]
		if !ok {
			return nil, fmt.Errorf("combine: missing parent model %q", sel.ParentName)
		}
		child, ok := byName[sel.ChildName]
		if !ok {
			return nil, fmt.Errorf("combine: missing child model %q", sel.ChildName)
		}

		newName := combine.CombinedModelName(parent.Name, child.Name)
		outPath := filepath.Join(workDir, fmt.Sprintf("combine_bitscore_round_%03d_%02d_%s.cm", round, i, combine.SafeFileStem(newName)))
		if _, err := combine.CombineParentWithTopChild(ctx, parent.Path, outPath, parent.Name, child.Name, realPool, realMatrix, topProportion); err != nil {
			return nil, fmt.Errorf("combine %s <- %s: %w", parent.Name, child.Name, err)
		}

		replaced[parent.Name] = pipeline.Model{Name: newName, Path: outPath}
		removed[child.Name] = true
	}

	next := make([]pipeline.Model, 0, len(models)-len(removed))
	for _, m := range models {
		if removed[m.Name] {
			continue
		}
		if repl, ok := replaced[m.Name]; ok {
			next = append(next, repl)
			continue
		}
		next = append(next, m)
	}
	return next, nil
}

// printOverview writes a summary of the run to w before any work starts:
// how many models, how many pairs will be reported (vs. the total
// possible), and the subprocess-call estimate from pipeline.Estimate --
// see that function's doc comment, and the top-level README's complexity
// discussion, for why call count (linear in model count) and actual
// scoring compute (quadratic) are reported/discussed separately rather
// than conflated into one number.
func printOverview(w *os.File, cfg pipeline.Config) {
	est := pipeline.Estimate(cfg)

	fmt.Fprintln(w, "=== cm-simile: bitscoresumdifference ===")
	fmt.Fprintf(w, "Models:           %d\n", est.NModels)
	names := make([]string, len(cfg.Models))
	for i, m := range cfg.Models {
		names[i] = m.Name
	}
	fmt.Fprintf(w, "  %s\n", strings.Join(names, ", "))
	if len(cfg.Query) > 0 {
		fmt.Fprintf(w, "Pairs to report:  %d (of %d total; query filter active)\n", est.NPairsReported, est.NPairsTotal)
	} else {
		fmt.Fprintf(w, "Pairs to report:  %d\n", est.NPairsTotal)
	}
	fmt.Fprintf(w, "npool=%d  nshuffles=%d (bootstrap resamples, cheap)  seed=%d  threads=%d\n",
		cfg.NPool, cfg.NShuffles, cfg.Seed, cfg.Threads)
	fmt.Fprintf(w, "reservoir: min=%d max=%d tolerance=%.4f\n",
		cfg.Reservoir.ReservoirMin, cfg.Reservoir.ReservoirMax, cfg.Reservoir.ReservoirTolerance)
	fmt.Fprintln(w, "")
	fmt.Fprintln(w, "Estimated subprocess calls (reservoir cost is a RANGE, not exact --")
	fmt.Fprintln(w, "growth is adaptive and stops once the convergence diagnostic")
	fmt.Fprintln(w, "stabilizes, which depends on the actual score data):")
	fmt.Fprintf(w, "  cmemit   (sampling):            %d\n", est.CmemitCalls)
	fmt.Fprintf(w, "  cmsearch (real pool):            %d\n", est.CmsearchCallsReal)
	fmt.Fprintf(w, "  cmsearch (null reservoir):        %d (best case, converges round 1) .. %d (worst case, hits reservoir-max)\n",
		est.CmsearchCallsReservoirBest, est.CmsearchCallsReservoirWorst)
	fmt.Fprintf(w, "  ------------------------------------------------\n")
	fmt.Fprintf(w, "  total cmsearch calls:            %d .. %d\n", est.CmsearchCallsTotalBest, est.CmsearchCallsTotalWorst)
	fmt.Fprintln(w, "")
	fmt.Fprintln(w, "Note: call count is linear in model count either way, but each call scores")
	fmt.Fprintln(w, "the whole pool for that stage, so total scoring work is O(models^2 x npool)")
	fmt.Fprintln(w, "-- nshuffles no longer multiplies this (bootstrap resampling from the")
	fmt.Fprintln(w, "reservoir is in-memory, not a subprocess call) -- see README for details.")
	fmt.Fprintln(w, "")
}

// progressBar renders a simple, dependency-free progress bar to stderr,
// one line per pipeline stage, overwritten in place via '\r'. Safe to
// call concurrently (pipeline.Config.OnProgress can be invoked from many
// goroutines at once, one per in-flight subprocess call).
type progressBar struct {
	w            *os.File
	mu           sync.Mutex
	currentStage string
}

var stageLabels = map[string]string{
	"sampling":     "Sampling sequences (cmemit)",
	"scoring_real": "Scoring real pool (cmsearch)",
}

// labelFor returns a human label for a pipeline stage. Reservoir growth
// stages are named "scoring_reservoir_round_N" (N not known upfront,
// since growth is adaptive -- see pipeline.Config.OnProgress's doc
// comment), so those get a prefix match rather than an exact-match
// lookup; anything else falls back to the raw stage string.
func labelFor(stage string) string {
	if label, ok := stageLabels[stage]; ok {
		return label
	}
	if strings.HasPrefix(stage, "scoring_reservoir_round_") {
		round := strings.TrimPrefix(stage, "scoring_reservoir_round_")
		return fmt.Sprintf("Scoring null reservoir, round %s (cmsearch)", round)
	}
	return stage
}

func newProgressBar(w *os.File) *progressBar {
	return &progressBar{w: w}
}

func (p *progressBar) onProgress(stage string, completed, total int) {
	p.mu.Lock()
	defer p.mu.Unlock()

	if p.currentStage != "" && p.currentStage != stage {
		fmt.Fprintln(p.w) // finish the previous stage's line
	}
	p.currentStage = stage

	label := labelFor(stage)

	const barWidth = 30
	frac := 0.0
	if total > 0 {
		frac = float64(completed) / float64(total)
	}
	filled := int(frac * barWidth)
	if filled > barWidth {
		filled = barWidth
	}
	bar := strings.Repeat("#", filled) + strings.Repeat("-", barWidth-filled)

	fmt.Fprintf(p.w, "\r%-36s [%s] %3d%% (%d/%d)", label, bar, int(frac*100), completed, total)
	if completed >= total {
		fmt.Fprintln(p.w)
	}
}

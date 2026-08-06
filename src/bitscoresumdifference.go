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
	"strings"
	"sync"

	"cmsimile/src/infernal"
	"cmsimile/src/output"
	"cmsimile/src/pipeline"
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
	nShuffles := flag.Int("nshuffles", 200, "dinucleotide-shuffle null replicates per pair")
	seed := flag.Int64("seed", 42, "base random seed (deterministic reproducibility across runs)")
	threads := flag.Int("threads", 4, "max concurrent cmemit/cmsearch subprocess calls dispatched at once")
	outPath := flag.String("out", "", "output TSV path (default: stdout)")

	cmemitBin := flag.String("cmemit-bin", "", "path to cmemit binary (default: resolve \"cmemit\" from PATH)")
	cmsearchBin := flag.String("cmsearch-bin", "", "path to cmsearch binary (default: resolve \"cmsearch\" from PATH)")
	cmstatBin := flag.String("cmstat-bin", "", "path to cmstat binary, used for -cm-file expansion (default: resolve \"cmstat\" from PATH)")
	cmfetchBin := flag.String("cmfetch-bin", "", "path to cmfetch binary, used for -cm-file expansion (default: resolve \"cmfetch\" from PATH)")

	cmsearchCPU := flag.Int("cmsearch-cpu", 1, "threads passed to cmsearch's own --cpu (internal to a single cmsearch call; multiplies with -threads)")
	minBitScore := flag.Float64("min-bitscore", -1000.0, "cmsearch -T inclusion threshold; must stay low enough that even weak/shuffled sequences are never excluded from scoring (see Scorer doc comment)")
	noFilters := flag.Bool("disable-filters", true, "pass --max to cmsearch, disabling acceleration filters so every sequence is fully scored by the CM (slower, but required for correct null-distribution statistics; see Scorer doc comment)")

	keepTemp := flag.Bool("keep-temp", false, "keep temporary FASTA/tblout files instead of deleting them (for debugging a parsing mismatch against real Infernal output)")
	quiet := flag.Bool("quiet", false, "suppress the run overview and progress bar (written to stderr; stdout/-out results are unaffected either way)")

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
		Seed:      *seed,
		Threads:   *threads,
		Emitter:   infernal.Emitter{BinPath: *cmemitBin},
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

	results, err := pipeline.Run(context.Background(), cfg)
	if err != nil {
		return err
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
	fmt.Fprintf(w, "npool=%d  nshuffles=%d  seed=%d  threads=%d\n",
		cfg.NPool, cfg.NShuffles, cfg.Seed, cfg.Threads)
	fmt.Fprintln(w, "")
	fmt.Fprintln(w, "Estimated subprocess calls:")
	fmt.Fprintf(w, "  cmemit   (sampling):        %d\n", est.CmemitCalls)
	fmt.Fprintf(w, "  cmsearch (real pool):       %d\n", est.CmsearchCallsReal)
	fmt.Fprintf(w, "  cmsearch (null replicates): %d  (= %d models x %d shuffles)\n", est.CmsearchCallsNull, est.NModels, cfg.NShuffles)
	fmt.Fprintf(w, "  ------------------------------------------------\n")
	fmt.Fprintf(w, "  total cmsearch calls:       %d\n", est.CmsearchCallsTotal)
	fmt.Fprintln(w, "")
	fmt.Fprintln(w, "Note: call count is linear in model count, but each call scores the")
	fmt.Fprintln(w, "whole pool for that stage, so total scoring work is O(models^2 x npool")
	fmt.Fprintln(w, "x (1+nshuffles)) -- see README for the full complexity discussion.")
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
	"scoring_null": "Scoring null replicates (cmsearch)",
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

	label := stageLabels[stage]
	if label == "" {
		label = stage
	}

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

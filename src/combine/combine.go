// Package combine holds the logic for combining covariance models by
// expanding one parent model with the top-scoring real sequences from a
// child model's pool.
package combine

import (
	"context"
	"fmt"
	"os"
	"path/filepath"
	"strings"

	"cmsimile/src/infernal"
	"cmsimile/src/matrix"
	"cmsimile/src/seqpool"
	"cmsimile/src/mafft"
)

// CombineResult reports what subset of the child model's real pool was
// used to expand the parent model.
type CombineResult struct {
	NChildTotal    int
	NChildSelected int
	SelectedIDs    []string
}

// CombinedModelName returns a deterministic merged model name.
func CombinedModelName(parent, child string) string {
	return parent + "__plus__" + child
}

// SafeFileStem makes a name safe for use as a filename stem.
func SafeFileStem(name string) string {
	replacer := strings.NewReplacer(" ", "_", "/", "_", "\\", "_", ":", "_", "\t", "_")
	return replacer.Replace(strings.TrimSpace(name))
}

// CombineParentWithTopChild expands parentModelPath by aligning only the
// top-scoring child real sequences (as scored against parentModelName)
// and rebuilding a new model at outputModelPath. cpu is passed straight
// through to calibration (cmcalibrate's own internal parallelism) --
// callers running several merges concurrently should divide their
// thread budget across calls rather than passing the full count to each.
func CombineParentWithTopChild(
	ctx context.Context,
	aModelPath string,
	outputModelPath string,
	aModelName string,
	bModelName string,
	realPool *seqpool.Pool,
	realScores *matrix.ScoreMatrix,
	topProportion float64,
	cpu int,
) (CombineResult, error) {
	if err := validateCombineInputs(realPool, realScores, topProportion, cpu); err != nil {
		return CombineResult{}, err
	}

	bSeqs := realPool.ByOrigin(bModelName)
	if len(bSeqs) == 0 {
		return CombineResult{}, fmt.Errorf("combine: no sequences found in real pool for model %q", bModelName)
	}
	aSeqs := realPool.ByOrigin(aModelName)
	if len(aSeqs) == 0 {
		return CombineResult{}, fmt.Errorf("combine: no sequences found in real pool for model %q", aModelName)
	}

	// top B sequences scored against model A
	topBIDs, err := realScores.GetTopScorers(topProportion, aModelName, seqIDs(bSeqs))
	if err != nil {
		return CombineResult{}, fmt.Errorf("combine: selecting top scorers: %w", err)
	}
	// top A sequences scored against model B
	topAIDs, err := realScores.GetTopScorers(topProportion, bModelName, seqIDs(aSeqs))
	if err != nil {
		return CombineResult{}, fmt.Errorf("combine: selecting top scorers: %w", err)
	}
	topIDs := append(topBIDs, topAIDs...)

	selectedSeqs := selectByID(topIDs, bSeqs, aSeqs)
	if len(selectedSeqs) == 0 {
		return CombineResult{}, fmt.Errorf("combine: top selection produced no sequences")
	}

	if err := os.MkdirAll(filepath.Dir(outputModelPath), 0o755); err != nil {
		return CombineResult{}, fmt.Errorf("combine: ensuring output directory: %w", err)
	}

	ctx = withCombineDefaults(ctx)

	alignmentPath, err := alignSelectedSequences(ctx, selectedSeqs)
	if err != nil {
		return CombineResult{}, err
	}
	defer os.Remove(alignmentPath)

	builder := infernal.Builder{}
	if err := builder.Build(ctx, alignmentPath, outputModelPath, false); err != nil {
		return CombineResult{}, fmt.Errorf("combine: rebuilding combined model: %w", err)
	}

	calibrator := infernal.Calibrator{CPU: cpu}
	if err := calibrator.Calibrate(ctx, outputModelPath); err != nil {
		return CombineResult{}, fmt.Errorf("combine: calibrating rebuilt model: %w", err)
	}

	return CombineResult{
		NChildTotal:    len(bSeqs),
		NChildSelected: len(selectedSeqs),
		SelectedIDs:    append([]string(nil), topIDs...),
	}, nil
}

func validateCombineInputs(pool *seqpool.Pool, scores *matrix.ScoreMatrix, topProportion float64, cpu int) error {
	if pool == nil {
		return fmt.Errorf("combine: realPool is nil")
	}
	if scores == nil {
		return fmt.Errorf("combine: realScores is nil")
	}
	if topProportion <= 0 || topProportion > 1 {
		return fmt.Errorf("combine: topProportion must be in (0,1], got %.4f", topProportion)
	}
	if cpu < 1 {
		return fmt.Errorf("combine: cpu must be >= 1, got %d", cpu)
	}
	return nil
}

func seqIDs(seqs []seqpool.Sequence) []string {
	ids := make([]string, len(seqs))
	for i, s := range seqs {
		ids[i] = s.ID
	}
	return ids
}

// selectByID returns, in the order the seqSets are given, every sequence
// whose ID appears in ids.
func selectByID(ids []string, seqSets ...[]seqpool.Sequence) []seqpool.Sequence {
	want := make(map[string]bool, len(ids))
	for _, id := range ids {
		want[id] = true
	}
	selected := make([]seqpool.Sequence, 0, len(ids))
	for _, seqs := range seqSets {
		for _, s := range seqs {
			if want[s.ID] {
				selected = append(selected, s)
			}
		}
	}
	return selected
}

// alignSelectedSequences writes seqs to a fresh temp alignment file and runs
// the aligner over it. On success, the caller owns the returned path and is
// responsible for removing it; on error, the temp file is cleaned up here.
func alignSelectedSequences(ctx context.Context, seqs []seqpool.Sequence) (path string, err error) {
	alignmentFile, err := os.CreateTemp("", "cmsimile-combine-*.sto")
	if err != nil {
		return "", fmt.Errorf("combine: creating temporary alignment path: %w", err)
	}
	alignmentPath := alignmentFile.Name()
	if err := alignmentFile.Close(); err != nil {
		os.Remove(alignmentPath)
		return "", fmt.Errorf("combine: closing temporary alignment file: %w", err)
	}

	aligner := mafft.Aligner{}
	if err := aligner.Align(ctx, seqs, alignmentPath); err != nil {
		os.Remove(alignmentPath)
		return "", fmt.Errorf("combine: aligning sequences: %w", err)
	}
	return alignmentPath, nil
}

func withCombineDefaults(ctx context.Context) context.Context {
	if ctx == nil {
		ctx = context.Background()
	}
	if _, ok := ctx.Value("seed").(int); !ok {
		ctx = context.WithValue(ctx, "seed", 1)
	}
	if _, ok := ctx.Value("cpu").(int); !ok {
		ctx = context.WithValue(ctx, "cpu", 1)
	}
	return ctx
}
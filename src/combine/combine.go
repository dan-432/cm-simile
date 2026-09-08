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
	parentModelPath string,
	outputModelPath string,
	parentModelName string,
	childModelName string,
	realPool *seqpool.Pool,
	realScores *matrix.ScoreMatrix,
	topProportion float64,
	cpu int,
) (CombineResult, error) {
	if realPool == nil {
		return CombineResult{}, fmt.Errorf("combine: realPool is nil")
	}
	if realScores == nil {
		return CombineResult{}, fmt.Errorf("combine: realScores is nil")
	}
	if topProportion <= 0 || topProportion > 1 {
		return CombineResult{}, fmt.Errorf("combine: topProportion must be in (0,1], got %.4f", topProportion)
	}
	if cpu < 1 {
		return CombineResult{}, fmt.Errorf("combine: cpu must be >= 1, got %d", cpu)
	}

	childSeqs := realPool.ByOrigin(childModelName)
	if len(childSeqs) == 0 {
		return CombineResult{}, fmt.Errorf("combine: no child sequences found in real pool for model %q", childModelName)
	}
	parentSeqs := realPool.ByOrigin(parentModelName)
	if len(parentSeqs) == 0 {
		return CombineResult{}, fmt.Errorf("combine: no parent sequences found in real pool for model %q", parentModelName)
	}

	alignIDs := make([]string, len(childSeqs)+len(parentSeqs))
	for i, s := range childSeqs {
		alignIDs[i] = s.ID
	}

	for i, s := range parentSeqs {
		alignIDs[len(childSeqs)+i] = s.ID
	}

	topIDs, err := realScores.GetTopScorers(topProportion, parentModelName, alignIDs)
	if err != nil {
		return CombineResult{}, fmt.Errorf("combine: selecting top scorers: %w", err)
	}

	selectedID := make(map[string]bool, len(topIDs))
	for _, id := range topIDs {
		selectedID[id] = true
	}

	selectedSeqs := make([]seqpool.Sequence, 0, len(topIDs))
	for _, s := range childSeqs {
		if selectedID[s.ID] {
			selectedSeqs = append(selectedSeqs, s)
		}
	}
	for _, s := range parentSeqs {
		if selectedID[s.ID] {
			selectedSeqs = append(selectedSeqs, s)
		}
	}
	if len(selectedSeqs) == 0 {
		return CombineResult{}, fmt.Errorf("combine: top child selection produced no sequences")
	}

	if err := os.MkdirAll(filepath.Dir(outputModelPath), 0o755); err != nil {
		return CombineResult{}, fmt.Errorf("combine: ensuring output directory: %w", err)
	}

	ctx = withCombineDefaults(ctx)
	alignmentFile, err := os.CreateTemp("", "cmsimile-combine-*.sto")
	if err != nil {
		return CombineResult{}, fmt.Errorf("combine: creating temporary alignment path: %w", err)
	}
	alignmentPath := alignmentFile.Name()
	alignmentFile.Close()
	defer os.Remove(alignmentPath)

	// aligner := infernal.Aligner{}
	aligner := mafft.Aligner{}
	if err := aligner.Align(ctx, selectedSeqs, alignmentPath); err != nil {
		return CombineResult{}, fmt.Errorf("combine: aligning child subset to parent: %w", err)
	}

	builder := infernal.Builder{}
	if err := builder.Build(ctx, alignmentPath, outputModelPath); err != nil {
		return CombineResult{}, fmt.Errorf("combine: rebuilding parent model: %w", err)
	}

	calibrator := infernal.Calibrator{CPU: cpu}
	if err := calibrator.Calibrate(ctx, outputModelPath); err != nil {
		return CombineResult{}, fmt.Errorf("combine: calibrating rebuilt model: %w", err)
	}

	return CombineResult{
		NChildTotal:    len(childSeqs),
		NChildSelected: len(selectedSeqs),
		SelectedIDs:    append([]string(nil), topIDs...),
	}, nil
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
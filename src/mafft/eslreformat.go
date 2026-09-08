package mafft

import (
	"context"
	"fmt"
)

// Reformatter wraps esl-reformat binary. It is a thin, mockable wrapper:
// production code uses Reformatter directly; tests can swap in a fake that
// satisfies the same method set (see ReformatterFunc-style fakes in pipeline
// tests TODO) without needing esl-reformat installed.
type Reformatter struct {
	// BinPath is the path to the esl-reformat binary. If empty, "esl-reformat" is
	// resolved from PATH at call time.
	BinPath string
}

func (r Reformatter) bin() string { return resolveBin(r.BinPath, "esl-reformat") }

func (r Reformatter) Reformat(ctx context.Context, inputPath, outputPath, format string) error {
	args := []string{
		format,
		inputPath,
	}

	err := runCapture(ctx, r.bin(), outputPath, args...)
	if err != nil {
		return fmt.Errorf("esl-reformat against %s: %w", inputPath, err)
	}
	return nil
}
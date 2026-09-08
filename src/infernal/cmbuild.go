package infernal

import (
	"context"
	"fmt"
)

type Builder struct {
	// BinPath is the path to the cmbuild binary. If empty, "cmbuild" is
	// resolved from PATH at call time.
	BinPath string
}

func (b Builder) bin() string { return resolveBin(b.BinPath, "cmbuild") }

func (b Builder) Build(ctx context.Context, alignmentPath, cmPath string) error {
	args := []string{
		"--hand",
		"--informat", "stockholm",
		cmPath,
		alignmentPath,
	}
	out, err := runCapture(ctx, b.bin(), args...)
	if err != nil {
		return fmt.Errorf("cmbuild from %s: %w\ncmbuild output:\n%s", alignmentPath, err, string(out))
	}
	return nil
}
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

func (b Builder) Build(ctx context.Context, alignmentPath, cmPath string, structured bool) error {
	// If structured is true, the --hand option is passed to cmbuild, which is for already structured alignments
	args := []string{}

	if structured {
		args = append(args, "--hand")
	} else {
		args = append(args, "--noss")
	}
	args = append(args, "--informat", "stockholm", cmPath, alignmentPath)

	out, err := runCapture(ctx, b.bin(), args...)
	if err != nil {
		return fmt.Errorf("cmbuild from %s: %w\ncmbuild output:\n%s", alignmentPath, err, string(out))
	}
	return nil
}
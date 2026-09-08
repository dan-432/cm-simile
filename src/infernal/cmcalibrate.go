package infernal

import (
	"context"
	"fmt"
	"strconv"
)

type Calibrator struct {
	// BinPath is the path to the cmcalibrate binary. If empty, "cmcalibrate" is
	// resolved from PATH at call time.
	BinPath string

	CPU int
}

func (c Calibrator) bin() string { return resolveBin(c.BinPath, "cmcalibrate") }

func (c Calibrator) Calibrate(ctx context.Context, cmPath string) error {
	args := []string{
		"--cpu", strconv.Itoa(c.CPU),
		"--seed", strconv.Itoa(ctx.Value("seed").(int)),
		cmPath,
	}
	out, err := runCapture(ctx, c.bin(), args...)
	if err != nil {
		return fmt.Errorf("cmcalibrate for %s: %w\ncmcalibrate output:\n%s", cmPath, err, string(out))
	}
	return nil
}
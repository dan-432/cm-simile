package infernal

import (
	"bytes"
	"context"
	"fmt"
	"os/exec"
)

// runCapture runs bin with args, returning stdout on success. stderr is
// included in the returned error on failure so a failing Infernal
// invocation (e.g. a malformed CM file) is diagnosable from the error
// message alone, without needing to re-run by hand.
func runCapture(ctx context.Context, bin string, args ...string) ([]byte, error) {
	cmd := exec.CommandContext(ctx, bin, args...)
	var stdout, stderr bytes.Buffer
	cmd.Stdout = &stdout
	cmd.Stderr = &stderr
	if err := cmd.Run(); err != nil {
		return nil, fmt.Errorf("running %s %v: %w\nstderr:\n%s", bin, args, err, stderr.String())
	}
	return stdout.Bytes(), nil
}

// resolveBin returns override if non-empty, otherwise the bare command
// name, letting exec.LookPath resolve it from PATH at run time.
func resolveBin(override, fallback string) string {
	if override != "" {
		return override
	}
	return fallback
}

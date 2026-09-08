package infernal

import (
	"bytes"
	"context"
	"fmt"
	"os"
	"strconv"

	"cmsimile/src/seqpool"
)

// Aligner wraps Infernal's cmalign binary. It is a thin, mockable wrapper:
// production code uses Aligner directly; tests can swap in a fake that
// satisfies the same method set (see AlignerFunc-style fakes in pipeline
// tests) without needing Infernal installed.
type Aligner struct {
	// BinPath is the path to the cmalign binary. If empty, "cmalign" is
	// resolved from PATH at call time.
	BinPath string
}

func (a Aligner) bin() string { return resolveBin(a.BinPath, "cmalign") }

func (a Aligner) Align(ctx context.Context, cmPath string, seqs []seqpool.Sequence, alignmentPath string) (err error) {
	if len(seqs) == 0 {
		return fmt.Errorf("infernal: Align called with empty sequence set")
	}

	tmpFasta, err := os.CreateTemp("", "cmsimile-pool-*.fasta")
	if err != nil {
		return fmt.Errorf("cmalign: creating temp FASTA: %w", err)
	}
	defer os.Remove(tmpFasta.Name())

	pool := seqpool.Pool{Sequences: seqs}
	var fastaBuf bytes.Buffer
	if err := pool.WriteFASTA(&fastaBuf); err != nil {
		return fmt.Errorf("cmalign: writing temp FASTA: %w", err)
	}
	_, err = tmpFasta.Write(fastaBuf.Bytes())
	if err != nil {
		return fmt.Errorf("cmalign: writing temp FASTA content: %w", err)
	}
	tmpFasta.Close()

	args := []string{
		"-o", alignmentPath,
		"--outformat", "Stockholm",
		"--cpu", strconv.Itoa(ctxIntValueOr(ctx, "cpu", 1)),
		cmPath,
		tmpFasta.Name(),
	}
	// Run cmalign and capture its output for debugging purposes
	out, err := runCapture(ctx, a.bin(), args...)
	if err != nil {
		return fmt.Errorf("cmalign against %s: %w\ncmalign output:\n%s", cmPath, err, string(out))
	}
	return nil
}

func ctxIntValueOr(ctx context.Context, key string, fallback int) int {
	v := ctx.Value(key)
	if v == nil {
		return fallback
	}
	iv, ok := v.(int)
	if !ok {
		return fallback
	}
	if iv <= 0 {
		return fallback
	}
	return iv
}

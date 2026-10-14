package mafft

import (
	// "bytes"
	"context"
	"fmt"
	"os"

	"cmsimile/src/seqpool"
)

// Aligner wraps mafft binary. It is a thin, mockable wrapper:
// production code uses Aligner directly; tests can swap in a fake that
// satisfies the same method set (see AlignerFunc-style fakes in pipeline
// tests TODO) without needing mafft installed.
type Aligner struct {
	// BinPath is the path to the mafft binary. If empty, "mafft" is
	// resolved from PATH at call time.
	BinPath string
}

func (a Aligner) bin() string { return resolveBin(a.BinPath, "mafft") }

func (a Aligner) Align(ctx context.Context, seqs []seqpool.Sequence, alignmentPath string) (err error) {
	if len(seqs) == 0 {
		return fmt.Errorf("mafft: Align called with empty sequence set")
	}

	tmpFasta, err := os.CreateTemp("", "cmsimile-pool-*.fasta")
	if err != nil {
		return fmt.Errorf("mafft: creating temp FASTA: %w", err)
	}
	defer os.Remove(tmpFasta.Name())

	// Write the sequences to the temporary FASTA file
	pool := seqpool.Pool{Sequences: seqs}

	if err := pool.WriteFASTA(tmpFasta); err != nil {
		return fmt.Errorf("mafft: writing temp FASTA: %w", err)
	}

	args := []string{
		"--localpair",
		"--maxiterate",
		"1000",
		tmpFasta.Name(),
	}

	tmpAlignFasta, err := os.CreateTemp("", "cmsimile-align-*.fasta")
	if err != nil {
		return fmt.Errorf("mafft: creating temp alignment FASTA: %w", err)
	}
	defer os.Remove(tmpAlignFasta.Name())
	if err := tmpAlignFasta.Close(); err != nil {
		return fmt.Errorf("mafft: closing temp alignment FASTA: %w", err)
	}

	// Run mafft and capture its output for debugging purposes, out is the stdout of mafft which is the aligned sequences in fasta format
	if err := runCapture(ctx, a.bin(), tmpAlignFasta.Name(), args...); err != nil {
		return fmt.Errorf("mafft against %s: %w", tmpFasta.Name(), err)
	}

	// reformat the aligned sequences to stockholm format using esl-reformat
	reformatter := Reformatter{}
	if err := reformatter.Reformat(ctx, tmpAlignFasta.Name(), alignmentPath, "stockholm"); err != nil {
		return fmt.Errorf("mafft: reformatting aligned sequences to stockholm format: %w", err)
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
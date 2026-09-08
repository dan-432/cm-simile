package mafft

import (
	"bytes"
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
	var fastaBuf bytes.Buffer
	if err := pool.WriteFASTA(&fastaBuf); err != nil {
		return fmt.Errorf("mafft: writing temp FASTA: %w", err)
	}
	_, err = tmpFasta.Write(fastaBuf.Bytes())
	if err != nil {
		return fmt.Errorf("mafft: writing temp FASTA content: %w", err)
	}
	tmpFasta.Close()

	args := []string{
		tmpFasta.Name(),
	}

	tmpAlignFasta, err := os.CreateTemp("", "cmsimile-align-*.fasta")
	if err != nil {
		return fmt.Errorf("mafft: creating temp alignment FASTA: %w", err)
	}
	defer os.Remove(tmpAlignFasta.Name())
	// Run mafft and capture its output for debugging purposes, out is the stdout of mafft which is the aligned sequences in fasta format
	err = runCapture(ctx, a.bin(), tmpAlignFasta.Name(), args...)
	if err != nil {
		return fmt.Errorf("mafft against %s: %w", tmpFasta.Name(), err)
	}

	// reformat the aligned sequences to stockholm format using esl-reformat
	reformatter := Reformatter{}
	err = reformatter.Reformat(ctx, tmpAlignFasta.Name(), alignmentPath, "stockholm")
	if err != nil {
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

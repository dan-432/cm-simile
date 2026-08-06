package infernal

import (
	"bytes"
	"context"
	"fmt"
	"strconv"

	"cmsimile/src/seqpool"
)

// Emitter draws sequences from a covariance model file via Infernal's
// cmemit. It is a thin, mockable wrapper: production code uses Emitter
// directly; tests can swap in a fake that satisfies the same method set
// (see EmitterFunc-style fakes in pipeline tests) without needing Infernal
// installed.
type Emitter struct {
	// BinPath is the path to the cmemit binary. If empty, "cmemit" is
	// resolved from PATH at call time.
	BinPath string
}

func (e Emitter) bin() string { return resolveBin(e.BinPath, "cmemit") }

// Sample draws n independent sequences from cmPath using cmemit's default
// (unaligned, generative) emission mode, seeded for reproducibility, and
// tags each with originLabel as its model of origin. IDs are generated as
// "<originLabel>_sample_NNNN" so they are unique within a pool built from
// multiple models.
//
// Command run: cmemit -N <n> --seed <seed> <cmPath>
func (e Emitter) Sample(ctx context.Context, cmPath, originLabel string, n int, seed int64) ([]seqpool.Sequence, error) {
	if n <= 0 {
		return nil, fmt.Errorf("infernal: Sample requires n > 0, got %d", n)
	}
	args := []string{"-N", strconv.Itoa(n), "--seed", strconv.FormatInt(seed, 10), cmPath}
	out, err := runCapture(ctx, e.bin(), args...)
	if err != nil {
		return nil, fmt.Errorf("cmemit sample from %s: %w", cmPath, err)
	}
	records, err := ParseFASTA(bytes.NewReader(out))
	if err != nil {
		return nil, fmt.Errorf("cmemit sample from %s: parsing output: %w", cmPath, err)
	}
	if len(records) != n {
		// cmemit can in principle emit fewer sequences than requested for
		// unusual/degenerate models. Surface this rather than silently
		// pooling an under-sized sample -- see the mean-vs-sum discussion:
		// downstream statistics assume matched pool sizes across models,
		// and a silent mismatch here is exactly the kind of thing that
		// should be a loud error, not a quiet bias.
		return nil, fmt.Errorf("cmemit sample from %s: requested %d sequences, got %d", cmPath, n, len(records))
	}
	seqs := make([]seqpool.Sequence, len(records))
	for i, rec := range records {
		id := fmt.Sprintf("%s_sample_%04d", originLabel, i)
		seqs[i] = seqpool.NewSequence(id, originLabel, rec.Seq)
	}
	return seqs, nil
}

// Consensus draws the single maximum-a-posteriori consensus sequence from
// cmPath, tagged as "<originLabel>_consensus".
//
// Command run: cmemit -c <cmPath>
func (e Emitter) Consensus(ctx context.Context, cmPath, originLabel string) (seqpool.Sequence, error) {
	out, err := runCapture(ctx, e.bin(), "-c", cmPath)
	if err != nil {
		return seqpool.Sequence{}, fmt.Errorf("cmemit consensus from %s: %w", cmPath, err)
	}
	records, err := ParseFASTA(bytes.NewReader(out))
	if err != nil {
		return seqpool.Sequence{}, fmt.Errorf("cmemit consensus from %s: parsing output: %w", cmPath, err)
	}
	if len(records) != 1 {
		return seqpool.Sequence{}, fmt.Errorf("cmemit consensus from %s: expected 1 record, got %d", cmPath, len(records))
	}
	id := originLabel + "_consensus"
	return seqpool.NewSequence(id, originLabel, records[0].Seq), nil
}

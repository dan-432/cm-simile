package stats

import (
	"testing"

	"cmsimile/src/matrix"
	"cmsimile/src/seqpool"
)

func TestComputePairStatHappyPath(t *testing.T) {
	pool := &seqpool.Pool{}
	pool.Add(
		seqpool.NewSequence("a1", "A", "AAAA"),
		seqpool.NewSequence("a2", "A", "CCCC"),
		seqpool.NewSequence("b1", "B", "GGGG"),
		seqpool.NewSequence("b2", "B", "TTTT"),
	)

	m := matrix.New()
	// A sequences
	m.Set("a1", "A", 10)
	m.Set("a2", "A", 20)
	m.Set("a1", "B", 8)
	m.Set("a2", "B", 12)
	// B sequences
	m.Set("b1", "B", 30)
	m.Set("b2", "B", 40)
	m.Set("b1", "A", 9)
	m.Set("b2", "A", 11)

	ps, err := ComputePairStat(pool, m, "A", "B", []float64{5, 8, 9}, []float64{9, 10, 12})
	if err != nil {
		t.Fatalf("ComputePairStat error: %v", err)
	}

	if ps.MedianScoreAonA != 15 || ps.MedianScoreAonB != 10 {
		t.Fatalf("unexpected A medians: AonA=%v AonB=%v", ps.MedianScoreAonA, ps.MedianScoreAonB)
	}
	if ps.MedianScoreBonB != 35 || ps.MedianScoreBonA != 10 {
		t.Fatalf("unexpected B medians: BonB=%v BonA=%v", ps.MedianScoreBonB, ps.MedianScoreBonA)
	}
	if ps.NShuffles != 3 {
		t.Fatalf("NShuffles=%d, want 3", ps.NShuffles)
	}
}

func TestComputePairStatMismatchedNullLengths(t *testing.T) {
	pool := &seqpool.Pool{}
	pool.Add(seqpool.NewSequence("a1", "A", "AAAA"), seqpool.NewSequence("b1", "B", "CCCC"))
	m := matrix.New()
	m.Set("a1", "A", 1)
	m.Set("a1", "B", 1)
	m.Set("b1", "A", 1)
	m.Set("b1", "B", 1)

	if _, err := ComputePairStat(pool, m, "A", "B", []float64{1, 2}, []float64{1}); err == nil {
		t.Fatalf("expected mismatch null length error")
	}
}

func TestEmpiricalPValue(t *testing.T) {
	if got := empiricalPValue(2.0, nil); got != 1.0 {
		t.Fatalf("empiricalPValue with empty null = %v, want 1", got)
	}
	got := empiricalPValue(3.0, []float64{1, 3, 5})
	// Two null values >= 3 (3 and 5), with +1/+1 correction => (2+1)/(3+1)=0.75
	if got != 0.75 {
		t.Fatalf("empiricalPValue = %v, want 0.75", got)
	}
}

package combine

import (
	"context"
	"testing"

	"cmsimile/src/matrix"
	"cmsimile/src/seqpool"
)

func TestCombinedModelNameAndSafeFileStem(t *testing.T) {
	if got := CombinedModelName("A", "B"); got != "A__plus__B" {
		t.Fatalf("CombinedModelName = %q, want A__plus__B", got)
	}
	if got := SafeFileStem(" A/B: C\\D\t"); got != "A_B__C_D" {
		t.Fatalf("SafeFileStem = %q, unexpected", got)
	}
}

func TestCombineParentWithTopChildValidation(t *testing.T) {
	pool := &seqpool.Pool{}
	scores := matrix.New()

	if _, err := CombineParentWithTopChild(context.Background(), "p.cm", "o.cm", "P", "C", nil, scores, 0.5, 4); err == nil {
		t.Fatalf("expected nil pool error")
	}
	if _, err := CombineParentWithTopChild(context.Background(), "p.cm", "o.cm", "P", "C", pool, nil, 0.5, 4); err == nil {
		t.Fatalf("expected nil score matrix error")
	}
	if _, err := CombineParentWithTopChild(context.Background(), "p.cm", "o.cm", "P", "C", pool, scores, 0, 4); err == nil {
		t.Fatalf("expected invalid proportion error")
	}
	if _, err := CombineParentWithTopChild(context.Background(), "p.cm", "o.cm", "P", "C", pool, scores, 1.5, 0); err == nil {
		t.Fatalf("expected invalid cpu error")
	}

	pool.Add(seqpool.NewSequence("c1", "C", "AAAA"))
	if _, err := CombineParentWithTopChild(context.Background(), "p.cm", "o.cm", "P", "C", pool, scores, 0.5, 4); err == nil {
		t.Fatalf("expected missing parent score error")
	}
}

func TestWithCombineDefaults(t *testing.T) {
	ctx := withCombineDefaults(context.Background())
	if got, ok := ctx.Value("seed").(int); !ok || got != 1 {
		t.Fatalf("seed default missing or wrong: %v (ok=%v)", got, ok)
	}
	if got, ok := ctx.Value("cpu").(int); !ok || got != 1 {
		t.Fatalf("cpu default missing or wrong: %v (ok=%v)", got, ok)
	}
}

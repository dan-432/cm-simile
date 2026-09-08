package matrix

import "testing"

func TestSetGetAndRowCopy(t *testing.T) {
	m := New()
	m.Set("s1", "mA", 1.5)
	m.Set("s1", "mB", 2.5)

	if got, ok := m.Get("s1", "mA"); !ok || got != 1.5 {
		t.Fatalf("Get(s1,mA) = (%v,%v), want (1.5,true)", got, ok)
	}

	row := m.Row("s1")
	if len(row) != 2 {
		t.Fatalf("Row(s1) len = %d, want 2", len(row))
	}
	row["mA"] = 999
	if got, _ := m.Get("s1", "mA"); got != 1.5 {
		t.Fatalf("Row returned internal map; got score %v after external mutation", got)
	}
}

func TestMeanMedianAndTopScorers(t *testing.T) {
	m := New()
	m.Set("a", "model", 10)
	m.Set("b", "model", 20)
	m.Set("c", "model", 30)
	m.Set("d", "model", 40)

	ids := []string{"a", "b", "c", "d"}

	mean, err := m.MeanFor(ids, "model")
	if err != nil {
		t.Fatalf("MeanFor returned error: %v", err)
	}
	if mean != 25 {
		t.Fatalf("MeanFor = %v, want 25", mean)
	}

	median, err := m.MedianFor(ids, "model")
	if err != nil {
		t.Fatalf("MedianFor returned error: %v", err)
	}
	if median != 25 {
		t.Fatalf("MedianFor = %v, want 25", median)
	}

	top, err := m.GetTopScorers(0.5, "model", ids)
	if err != nil {
		t.Fatalf("GetTopScorers returned error: %v", err)
	}
	if len(top) != 2 {
		t.Fatalf("GetTopScorers len = %d, want 2", len(top))
	}
	if top[0] != "d" || top[1] != "c" {
		t.Fatalf("GetTopScorers order = %v, want [d c]", top)
	}
}

func TestErrorsForMissingAndInvalidInput(t *testing.T) {
	m := New()
	m.Set("a", "model", 1)

	if _, err := m.MustGet("missing", "model"); err == nil {
		t.Fatalf("MustGet on missing entry expected error")
	}
	if _, err := m.MeanFor(nil, "model"); err == nil {
		t.Fatalf("MeanFor(nil) expected error")
	}
	if _, err := m.MedianFor(nil, "model"); err == nil {
		t.Fatalf("MedianFor(nil) expected error")
	}
	if _, err := m.GetTopScorers(0, "model", []string{"a"}); err == nil {
		t.Fatalf("GetTopScorers proportion 0 expected error")
	}
	if _, err := m.GetTopScorers(1.1, "model", []string{"a"}); err == nil {
		t.Fatalf("GetTopScorers proportion >1 expected error")
	}
}

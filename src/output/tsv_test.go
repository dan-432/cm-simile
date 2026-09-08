package output

import (
	"bytes"
	"strings"
	"testing"

	"cmsimile/src/stats"
)

func TestWriteTSV_HeaderAndRowFormatting(t *testing.T) {
	rows := []stats.PairStat{{
		ModelA: "A", ModelB: "B",
		MedianScoreAonA: 1.23456, MedianScoreAonB: 2.34567,
		MedianScoreBonB: 3.45678, MedianScoreBonA: 4.56789,
		StatAintoB: 5.67891, StatBintoA: 6.78912,
		PValueAintoB: 0.01234, PValueBintoA: 0.05678,
		NShuffles: 42,
	}}

	var buf bytes.Buffer
	if err := WriteTSV(&buf, rows); err != nil {
		t.Fatalf("WriteTSV error: %v", err)
	}

	out := buf.String()
	lines := strings.Split(strings.TrimSpace(out), "\n")
	if len(lines) != 2 {
		t.Fatalf("expected 2 lines, got %d", len(lines))
	}
	if !strings.Contains(lines[0], "model_a\tmodel_b") {
		t.Fatalf("missing expected header, got: %s", lines[0])
	}
	if !strings.HasPrefix(lines[1], "A\tB\t") {
		t.Fatalf("data row prefix mismatch: %s", lines[1])
	}
	if !strings.Contains(lines[1], "\t42") {
		t.Fatalf("n_shuffles missing from row: %s", lines[1])
	}
}

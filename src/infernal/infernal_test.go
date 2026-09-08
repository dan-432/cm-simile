package infernal

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
)

func TestParseFASTA(t *testing.T) {
	text := ">h1\nAC\nTG\n\n>h2\nUUUU\n"
	recs, err := ParseFASTA(strings.NewReader(text))
	if err != nil {
		t.Fatalf("ParseFASTA error: %v", err)
	}
	if len(recs) != 2 {
		t.Fatalf("records len = %d, want 2", len(recs))
	}
	if recs[0].Header != "h1" || recs[0].Seq != "ACTG" {
		t.Fatalf("record 0 = %+v, want Header=h1 Seq=ACTG", recs[0])
	}
	if recs[1].Header != "h2" || recs[1].Seq != "UUUU" {
		t.Fatalf("record 1 = %+v, want Header=h2 Seq=UUUU", recs[1])
	}
}

func TestParseFASTARejectsDataBeforeHeader(t *testing.T) {
	if _, err := ParseFASTA(strings.NewReader("ACGT\n>h\nTT\n")); err == nil {
		t.Fatalf("expected error for sequence data before header")
	}
}

func TestParseTbloutBestHitPerSequence(t *testing.T) {
	tmp, err := os.CreateTemp("", "tblout-*.txt")
	if err != nil {
		t.Fatalf("CreateTemp: %v", err)
	}
	defer os.Remove(tmp.Name())

	content := "# target name         accession query name           accession mdl mdl from   mdl to       seq from         seq to strand trunc pass   gc  bias  score   E-value inc description of target\n"+
		"# ------------------- --------- -------------------- --------- --- -------- -------- -------------- -------------- ------ ----- ---- ---- ----- ------ --------- --- ---------------------\n"+
		"seqA     -         5S_rRNA              RF00001    cm        1      119              1            121      +    no    1 0.61   0.0  108.2   1.5e-27 ! -\n"+
		"seqA      -         tRNA                 RF00005    cm        1       73        2451092        2451021      -    no    1 0.55   0.0   64.3   2.1e-09 ! Transfer RNA tyrosine\n"+
		"seqB         -         U2                   RF00004    cm       10      150            45            189      +   5'    2 0.42   0.2   45.1     0.004 ? Truncated variant\n"+
		"# [ok]"
	if _, err := tmp.WriteString(content); err != nil {
		t.Fatalf("write temp tblout: %v", err)
	}
	if err := tmp.Close(); err != nil {
		t.Fatalf("close temp tblout: %v", err)
	}

	scores, err := parseTblout(tmp.Name())
	if err != nil {
		t.Fatalf("parseTblout error: %v", err)
	}
	if scores["seqA"] != 108.2 {
		t.Fatalf("seqA best score = %v, want 108.2", scores["seqA"])
	}
	if scores["seqB"] != 45.1 {
		t.Fatalf("seqB score = %v, want 45.1", scores["seqB"])
	}
}

func TestCopyFile(t *testing.T) {
	dir := t.TempDir()
	src := filepath.Join(dir, "in.txt")
	dst := filepath.Join(dir, "out.txt")
	if err := os.WriteFile(src, []byte("hello"), 0o644); err != nil {
		t.Fatalf("write src: %v", err)
	}
	if err := copyFile(src, dst); err != nil {
		t.Fatalf("copyFile: %v", err)
	}
	b, err := os.ReadFile(dst)
	if err != nil {
		t.Fatalf("read dst: %v", err)
	}
	if string(b) != "hello" {
		t.Fatalf("copied content = %q, want hello", string(b))
	}
}

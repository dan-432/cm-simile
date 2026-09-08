package seqpool

import (
	"bytes"
	"strings"
	"testing"
)

func TestNewSequenceNormalizesRNA(t *testing.T) {
	s := NewSequence("id", "orig", "augu")
	if s.Residues != "ATGT" {
		t.Fatalf("Residues = %q, want ATGT", s.Residues)
	}
}

func TestPoolAddGetByOriginAndOrigins(t *testing.T) {
	p := &Pool{}
	p.Add(
		NewSequence("a", "m1", "AAAA"),
		NewSequence("b", "m1", "CCCC"),
		NewSequence("c", "m2", "GGGG"),
	)

	if p.Len() != 3 {
		t.Fatalf("Len = %d, want 3", p.Len())
	}

	m1 := p.ByOrigin("m1")
	if len(m1) != 2 {
		t.Fatalf("ByOrigin(m1) len = %d, want 2", len(m1))
	}

	origins := p.Origins()
	if len(origins) != 2 || origins[0] != "m1" || origins[1] != "m2" {
		t.Fatalf("Origins = %v, want [m1 m2]", origins)
	}

	if seq, ok := p.Get("b"); !ok || seq.Origin != "m1" {
		t.Fatalf("Get(b) = (%v,%v), want found origin m1", seq, ok)
	}
	if _, ok := p.Get("missing"); ok {
		t.Fatalf("Get(missing) expected ok=false")
	}
}

func TestWriteFASTAAndDuplicateIDError(t *testing.T) {
	p := &Pool{}
	p.Add(NewSequence("a", "m", "ACTG"), NewSequence("b", "m", "TTTT"))

	var buf bytes.Buffer
	if err := p.WriteFASTA(&buf); err != nil {
		t.Fatalf("WriteFASTA error: %v", err)
	}
	out := buf.String()
	if !strings.Contains(out, ">a\nACTG\n") || !strings.Contains(out, ">b\nTTTT\n") {
		t.Fatalf("unexpected FASTA output:\n%s", out)
	}

	pDup := &Pool{}
	pDup.Add(NewSequence("dup", "m1", "AAAA"), NewSequence("dup", "m2", "CCCC"))
	if err := pDup.WriteFASTA(&bytes.Buffer{}); err == nil {
		t.Fatalf("expected duplicate ID error")
	}
}

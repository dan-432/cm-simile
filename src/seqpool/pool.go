// Package seqpool defines the core data structure cm-simile is built
// around: a flat, tagged pool of sequences drawn from one or more
// covariance models.
//
// The pool/matrix formulation (as opposed to a pairwise-loop formulation)
// is what lets cm-simile generalise from 1-to-1 comparison to bipartite
// (query set vs. reference set) and full many-to-many comparison without
// changing the scoring logic: every sequence, regardless of which model it
// came from, lives in one pool; every model scores that whole pool in a
// single cmscore call; pairwise statistics are then computed by slicing
// the resulting (sequence x model) score matrix in memory. See
// internal/matrix and internal/pipeline for how the pool is consumed.
package seqpool

import (
	"fmt"
	"io"
	"strings"
)

// Sequence is a single nucleotide sequence, tagged with the covariance
// model it was generated from. Origin is what lets many models' sequences
// share one pool without losing track of provenance.
type Sequence struct {
	// ID is a unique identifier within a Pool, e.g. "motif17_sample_0003"
	// or "motif17_consensus". Used as the FASTA header and as the key
	// into ScoreMatrix rows.
	ID string

	// Origin is the name of the covariance model this sequence was
	// sampled or derived from. Multiple sequences share the same Origin.
	Origin string

	// Residues is the nucleotide sequence itself: uppercase, ungapped,
	// alphabet {A,C,G,T} (U is normalised to T on construction so that
	// downstream dinucleotide shuffling and scoring don't need to special
	// case RNA vs DNA alphabets).
	Residues string
}

// NewSequence constructs a Sequence, normalising residues to uppercase DNA
// alphabet (T rather than U) so every downstream consumer can assume a
// consistent representation.
func NewSequence(id, origin, residues string) Sequence {
	r := strings.ToUpper(residues)
	r = strings.ReplaceAll(r, "U", "T")
	return Sequence{ID: id, Origin: origin, Residues: r}
}

// Pool is a flat collection of tagged sequences, potentially drawn from
// many different covariance models. It is the unit that gets written to
// FASTA and handed to cmscore: one pool, scored once per model, rather
// than one FASTA file per pair.
type Pool struct {
	Sequences []Sequence

	// index supports O(1) lookup by ID; built lazily by ensureIndex.
	index map[string]int
}

// Add appends sequences to the pool and invalidates the ID index.
func (p *Pool) Add(seqs ...Sequence) {
	p.Sequences = append(p.Sequences, seqs...)
	p.index = nil
}

// ByOrigin returns all sequences in the pool tagged with the given origin
// model name, in the order they were added.
func (p *Pool) ByOrigin(origin string) []Sequence {
	out := make([]Sequence, 0, len(p.Sequences))
	for _, s := range p.Sequences {
		if s.Origin == origin {
			out = append(out, s)
		}
	}
	return out
}

// Origins returns the distinct origin model names present in the pool, in
// first-seen order.
func (p *Pool) Origins() []string {
	seen := make(map[string]bool)
	var out []string
	for _, s := range p.Sequences {
		if !seen[s.Origin] {
			seen[s.Origin] = true
			out = append(out, s.Origin)
		}
	}
	return out
}

// Get looks up a sequence by ID. Reports ok=false if not present.
func (p *Pool) Get(id string) (Sequence, bool) {
	p.ensureIndex()
	i, ok := p.index[id]
	if !ok {
		return Sequence{}, false
	}
	return p.Sequences[i], true
}

func (p *Pool) ensureIndex() {
	if p.index != nil {
		return
	}
	p.index = make(map[string]int, len(p.Sequences))
	for i, s := range p.Sequences {
		p.index[s.ID] = i
	}
}

// WriteFASTA writes every sequence in the pool as a single multi-FASTA
// stream, suitable as direct input to cmscore. IDs must be unique within
// the pool; WriteFASTA returns an error if a duplicate is found so a
// silent ID collision can never quietly corrupt downstream score-matrix
// lookups.
func (p *Pool) WriteFASTA(w io.Writer) error {
	seen := make(map[string]bool, len(p.Sequences))
	for _, s := range p.Sequences {
		if seen[s.ID] {
			return fmt.Errorf("seqpool: duplicate sequence ID %q would corrupt FASTA/score-matrix correspondence", s.ID)
		}
		seen[s.ID] = true
		if _, err := fmt.Fprintf(w, ">%s\n%s\n", s.ID, s.Residues); err != nil {
			return err
		}
	}
	return nil
}

// Len reports the number of sequences in the pool.
func (p *Pool) Len() int { return len(p.Sequences) }

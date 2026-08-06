// Package shuffle implements the Altschul-Erikson dinucleotide-preserving
// sequence shuffle, used to build the empirical null distributions
// cm-simile's p-values are computed against.
package shuffle

import (
	"fmt"
	"math/rand"
	"strings"
)

// DinucleotideShuffle returns a randomised permutation of seq that
// preserves the first symbol, the last symbol, and the exact
// dinucleotide (2-mer) frequency distribution of the input, while
// destroying higher-order structure -- base-pairing complementarity,
// motif patterns, anything beyond local doublet composition.
//
// This is the null model cm-simile bootstraps against a covariance
// model's bit score on a dinucleotide-shuffled sequence reflects what the
// model would score for "typical composition, no real structural signal",
// isolating genuine shared sequence/structure signal between two models
// from a shared compositional bias they might independently favour (e.g.
// both preferring GC-rich sequence).
//
// Algorithm (Altschul & Erikson 1985; see also Clote et al. 2005 for the
// standard presentation): build the directed multigraph of observed
// dinucleotide transitions. For each symbol, fix its LAST outgoing
// transition (in the order it occurs in seq) as a distinguished edge. This
// specific selection -- the last exit from each symbol on a genuine walk
// through the sequence -- is guaranteed to form a tree of edges oriented
// towards the sequence's final symbol, because it is drawn from an actual
// Eulerian walk over that same graph. The remaining (non-fixed) outgoing
// edges of each symbol are randomly permuted, the fixed edge is appended
// last, and a new sequence is reconstructed by walking the shuffled edge
// lists starting from the original first symbol. The result is drawn
// uniformly from the set of all sequences sharing the same dinucleotide
// content as the input.
//
// seq is upper-cased and has U normalised to T before processing (see
// seqpool.NewSequence, which callers are expected to have already applied;
// this function re-applies the normalisation defensively so it behaves
// correctly even if called directly). Sequences shorter than 2 symbols are
// returned unchanged, since there is no dinucleotide content to preserve.
func DinucleotideShuffle(seq string, rng *rand.Rand) (string, error) {
	seq = strings.ToUpper(seq)
	seq = strings.ReplaceAll(seq, "U", "T")
	n := len(seq)
	if n < 2 {
		return seq, nil
	}

	// successors[x] = ordered list of symbols observed immediately after
	// x in seq, in first-to-last occurrence order.
	successors := make(map[byte][]byte)
	for i := 0; i < n-1; i++ {
		x, y := seq[i], seq[i+1]
		successors[x] = append(successors[x], y)
	}

	last := seq[n-1]

	// fixedEdge[x]: the LAST successor observed for x in the original
	// sequence. This becomes x's final outgoing edge after shuffling.
	fixedEdge := make(map[byte]byte, len(successors))
	for x, ys := range successors {
		fixedEdge[x] = ys[len(ys)-1]
		successors[x] = ys[:len(ys)-1] // remove it from the shuffle pool
	}

	// Defensive invariant check. Given fixedEdge is derived from a genuine
	// walk through seq, every symbol that appears as a source must be
	// able to reach `last` by following fixed edges. This should always
	// hold by construction; treat a violation as a programming error
	// (e.g. a bug in the code above) rather than proceeding to produce a
	// silently-wrong shuffle.
	if err := checkConnectedToLast(fixedEdge, last); err != nil {
		return "", fmt.Errorf("shuffle: internal invariant violated: %w", err)
	}

	// Randomly permute each symbol's remaining out-edges (Fisher-Yates),
	// then re-append its fixed last edge.
	for x, ys := range successors {
		rng.Shuffle(len(ys), func(i, j int) { ys[i], ys[j] = ys[j], ys[i] })
		successors[x] = append(ys, fixedEdge[x])
	}

	// Reconstruct the sequence: walk from the original first symbol,
	// consuming one successor at a time from each symbol's shuffled list.
	cursor := make(map[byte]int, len(successors))
	var b strings.Builder
	b.Grow(n)
	b.WriteByte(seq[0])
	cur := seq[0]
	for i := 0; i < n-1; i++ {
		idx := cursor[cur]
		list := successors[cur]
		if idx >= len(list) {
			return "", fmt.Errorf("shuffle: internal error: exhausted successor list for %q at position %d", string(cur), i)
		}
		next := list[idx]
		cursor[cur] = idx + 1
		b.WriteByte(next)
		cur = next
	}
	return b.String(), nil
}

// checkConnectedToLast verifies that every symbol appearing as a source in
// fixedEdge can reach `last` by following fixed edges (i.e. the fixed
// edges form a tree oriented towards the terminal symbol). This is an
// invariant that should always hold for fixedEdge derived from a genuine
// walk; see the algorithm note in DinucleotideShuffle's doc comment.
func checkConnectedToLast(fixedEdge map[byte]byte, last byte) error {
	reachable := map[byte]bool{last: true}
	for changed := true; changed; {
		changed = false
		for x, y := range fixedEdge {
			if reachable[y] && !reachable[x] {
				reachable[x] = true
				changed = true
			}
		}
	}
	for x := range fixedEdge {
		if x != last && !reachable[x] {
			return fmt.Errorf("symbol %q cannot reach terminal symbol %q via fixed edges", string(x), string(last))
		}
	}
	return nil
}

package shuffle

import (
	"math/rand"
	"testing"
)

// dinucCounts returns a map of every 2-mer's frequency in seq.
func dinucCounts(seq string) map[string]int {
	counts := make(map[string]int)
	for i := 0; i < len(seq)-1; i++ {
		counts[seq[i:i+2]]++
	}
	return counts
}

func mapsEqual(a, b map[string]int) bool {
	if len(a) != len(b) {
		return false
	}
	for k, v := range a {
		if b[k] != v {
			return false
		}
	}
	return true
}

func TestDinucleotideShuffle_PreservesInvariants(t *testing.T) {
	seqs := []string{
		"ACGUACGUACGUACGUACGU",
		"AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT",
		"GCGCGCGCGCGCGCGCGCGCATATATATATATATATATAT",
		"ATCGATCGATCGATCGGGGGCCCCCAAAAATTTTTACGT",
	}

	rng := rand.New(rand.NewSource(42))

	for _, seq := range seqs {
		// normalise the same way DinucleotideShuffle does internally, so
		// the "before" counts are computed on the same alphabet as the
		// "after" counts.
		normalized := normalizeForTest(seq)
		wantCounts := dinucCounts(normalized)
		wantFirst, wantLast := normalized[0], normalized[len(normalized)-1]

		for trial := 0; trial < 20; trial++ {
			shuffled, err := DinucleotideShuffle(seq, rng)
			if err != nil {
				t.Fatalf("seq %q trial %d: unexpected error: %v", seq, trial, err)
			}
			if len(shuffled) != len(normalized) {
				t.Fatalf("seq %q trial %d: length changed: got %d want %d", seq, trial, len(shuffled), len(normalized))
			}
			if shuffled[0] != wantFirst {
				t.Errorf("seq %q trial %d: first symbol changed: got %q want %q", seq, trial, string(shuffled[0]), string(wantFirst))
			}
			if shuffled[len(shuffled)-1] != wantLast {
				t.Errorf("seq %q trial %d: last symbol changed: got %q want %q", seq, trial, string(shuffled[len(shuffled)-1]), string(wantLast))
			}
			gotCounts := dinucCounts(shuffled)
			if !mapsEqual(wantCounts, gotCounts) {
				t.Errorf("seq %q trial %d: dinucleotide counts changed:\n before: %v\n after:  %v", seq, trial, wantCounts, gotCounts)
			}
		}
	}
}

func TestDinucleotideShuffle_ActuallyShuffles(t *testing.T) {
	// A sanity check that the output isn't trivially just the input
	// unchanged (which would technically satisfy the invariant tests
	// above but defeat the purpose). Use a sequence long/varied enough
	// that an identity "shuffle" across many trials would be extremely
	// unlikely by chance.
	seq := "ACGTACGTGGCATCGATCGATTAGCGATCGATCGGGCTATCGATCGATCG"
	rng := rand.New(rand.NewSource(7))

	differed := false
	for trial := 0; trial < 20; trial++ {
		shuffled, err := DinucleotideShuffle(seq, rng)
		if err != nil {
			t.Fatalf("unexpected error: %v", err)
		}
		if shuffled != normalizeForTest(seq) {
			differed = true
			break
		}
	}
	if !differed {
		t.Error("shuffle output was identical to input across 20 trials; expected genuine permutation")
	}
}

func TestDinucleotideShuffle_ShortSequenceUnchanged(t *testing.T) {
	rng := rand.New(rand.NewSource(1))
	for _, seq := range []string{"", "A", "AC"} {
		got, err := DinucleotideShuffle(seq, rng)
		if err != nil {
			t.Fatalf("seq %q: unexpected error: %v", seq, err)
		}
		want := normalizeForTest(seq)
		if len(seq) < 2 {
			if got != want {
				t.Errorf("seq %q: got %q want %q (should be unchanged, only normalised)", seq, got, want)
			}
		}
	}
}

func normalizeForTest(seq string) string {
	out := make([]byte, len(seq))
	for i := 0; i < len(seq); i++ {
		c := seq[i]
		if c >= 'a' && c <= 'z' {
			c -= 'a' - 'A'
		}
		if c == 'U' {
			c = 'T'
		}
		out[i] = c
	}
	return string(out)
}

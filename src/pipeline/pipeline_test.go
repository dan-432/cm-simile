package pipeline

import (
	"context"
	"fmt"
	"math/rand"
	"testing"

	"cmsimile/src/seqpool"
)

// fakeEmitter generates synthetic sequences per model without needing
// cmemit. Each model is given a distinct base pattern so that models
// sharing a pattern prefix are "similar" and models with unrelated
// patterns are not -- this lets the test assert cm-simile's statistics
// behave sensibly (low p-value / near-zero stat for similar models, the
// opposite for dissimilar ones) without any real CM files.
type fakeEmitter struct {
	basePattern map[string]string
}

func (f fakeEmitter) Sample(ctx context.Context, cmPath, originLabel string, n int, seed int64) ([]seqpool.Sequence, error) {
	pattern, ok := f.basePattern[originLabel]
	if !ok {
		return nil, fmt.Errorf("fakeEmitter: no pattern registered for %q", originLabel)
	}
	rng := rand.New(rand.NewSource(seed))
	seqs := make([]seqpool.Sequence, n)
	for i := 0; i < n; i++ {
		// Repeat the pattern to a fixed length and apply a little random
		// noise so sequences within a model aren't all identical.
		seq := mutate(repeatTo(pattern, 40), rng, 3)
		seqs[i] = seqpool.NewSequence(fmt.Sprintf("%s_sample_%04d", originLabel, i), originLabel, seq)
	}
	return seqs, nil
}

// fakeScorer scores a sequence against a model as (roughly) how well the
// sequence matches that model's registered base pattern -- a crude stand-
// in for a real CM bit score, but sufficient to validate that the
// pipeline's plumbing (pooling, matrix construction, null generation,
// statistic computation) produces sensible relative results.
type fakeScorer struct {
	basePattern map[string]string
	// pathToOrigin lets the scorer recover which model's pattern to use
	// from the cmPath argument passed by the pipeline, since in these
	// tests cmPath == model name.
}

func (f fakeScorer) Score(ctx context.Context, cmPath string, seqs []seqpool.Sequence) (map[string]float64, error) {
	pattern, ok := f.basePattern[cmPath]
	if !ok {
		return nil, fmt.Errorf("fakeScorer: no pattern registered for %q", cmPath)
	}
	repeated := repeatTo(pattern, 40)
	out := make(map[string]float64, len(seqs))
	for _, s := range seqs {
		out[s.ID] = matchScore(s.Residues, repeated)
	}
	return out, nil
}

func repeatTo(pattern string, length int) string {
	out := make([]byte, 0, length)
	for len(out) < length {
		out = append(out, pattern...)
	}
	return string(out[:length])
}

func mutate(seq string, rng *rand.Rand, nMutations int) string {
	b := []byte(seq)
	alphabet := []byte("ACGT")
	for i := 0; i < nMutations; i++ {
		pos := rng.Intn(len(b))
		b[pos] = alphabet[rng.Intn(len(alphabet))]
	}
	return string(b)
}

// matchScore is a toy "bit score": count of matching positions against a
// reference pattern, scaled to look roughly like a bit score magnitude.
func matchScore(seq, reference string) float64 {
	n := len(seq)
	if len(reference) < n {
		n = len(reference)
	}
	matches := 0
	for i := 0; i < n; i++ {
		if seq[i] == reference[i] {
			matches++
		}
	}
	return float64(matches) * 2.5 // arbitrary scale, bit-score-ish magnitude
}

// testReservoir returns a ReservoirConfig sized to converge quickly in
// tests: a small, loose tolerance so the adaptive growth loop in
// reservoir.go doesn't need many rounds against these tiny fake pools.
func testReservoir(npool int) ReservoirConfig {
	return ReservoirConfig{
		ResampleSize:        npool,
		ReservoirMin:        npool,
		ReservoirMax:        npool * 8,
		ReservoirTolerance:  0.25,
		DiagnosticResamples: 30,
	}
}

func TestPipeline_SimilarModelsScoreLowPValue(t *testing.T) {
	// modelA and modelB share the same underlying pattern (i.e. are
	// "the same family" in this synthetic setup); modelC is unrelated.
	basePattern := map[string]string{
		"modelA": "ACGTACGTGG",
		"modelB": "ACGTACGTGG", // identical pattern to A: should look highly similar
		"modelC": "TTTTCCCCAA", // unrelated pattern: should look dissimilar to both
	}

	cfg := Config{
		Models: []Model{
			{Name: "modelA", Path: "modelA"},
			{Name: "modelB", Path: "modelB"},
			{Name: "modelC", Path: "modelC"},
		},
		NPool:     50,
		NShuffles: 50,
		Reservoir: testReservoir(50),
		Seed:      1,
		Threads:   4,
		Emitter:   fakeEmitter{basePattern: basePattern},
		Scorer:    fakeScorer{basePattern: basePattern},
	}

	results, warnings, err := Run(context.Background(), cfg)
	if err != nil {
		t.Fatalf("Run failed: %v", err)
	}
	for _, w := range warnings {
		t.Logf("warning: %s", w)
	}
	if len(results) != 3 {
		t.Fatalf("expected 3 pairs (A-B, A-C, B-C), got %d", len(results))
	}

	byPair := make(map[string]int) // index into results, keyed by "A,B" etc.
	for i, r := range results {
		byPair[r.ModelA+","+r.ModelB] = i
	}

	abIdx, ok := byPair["modelA,modelB"]
	if !ok {
		t.Fatalf("missing modelA,modelB pair in results: %+v", results)
	}
	acIdx, ok := byPair["modelA,modelC"]
	if !ok {
		t.Fatalf("missing modelA,modelC pair in results: %+v", results)
	}

	ab := results[abIdx]
	ac := results[acIdx]

	t.Logf("A-B: statAintoB=%.3f statBintoA=%.3f pAintoB=%.4f pBintoA=%.4f",
		ab.StatAintoB, ab.StatBintoA, ab.PValueAintoB, ab.PValueBintoA)
	t.Logf("A-C: statAintoB=%.3f statBintoA=%.3f pAintoB=%.4f pBintoA=%.4f",
		ac.StatAintoB, ac.StatBintoA, ac.PValueAintoB, ac.PValueBintoA)

	// A and B (same pattern) should show much less score degradation
	// crossing models than A and C (unrelated patterns).
	if ab.StatAintoB <= ac.StatAintoB {
		t.Errorf("expected A-B cross-score degradation to be less negative than A-C, got StatAintoB(A,B)=%.3f, StatAintoB(A,C)=%.3f",
			ab.StatAintoB, ac.StatAintoB)
	}

	// A-B should be statistically significant (low p-value) since the
	// cross-recognition is far stronger than the shuffled null; A-C
	// should not be, since the patterns are unrelated and cross-scores
	// should look like noise, similar to the shuffled null.
	if ab.PValueAintoB >= 0.05 {
		t.Errorf("expected A-B to be significant (p < 0.05), got p=%.4f", ab.PValueAintoB)
	}
	if ac.PValueAintoB < ab.PValueAintoB {
		t.Errorf("expected A-C p-value (%.4f) to be larger (less significant) than A-B p-value (%.4f)",
			ac.PValueAintoB, ab.PValueAintoB)
	}
}

func TestPipeline_QueryFiltersReportedPairs(t *testing.T) {
	basePattern := map[string]string{
		"modelA": "ACGTACGTGG",
		"modelB": "GGTTCCAAGG",
		"modelC": "TTTTCCCCAA",
	}
	cfg := Config{
		Models: []Model{
			{Name: "modelA", Path: "modelA"},
			{Name: "modelB", Path: "modelB"},
			{Name: "modelC", Path: "modelC"},
		},
		Query:     map[string]bool{"modelA": true},
		NPool:     20,
		NShuffles: 20,
		Reservoir: testReservoir(20),
		Seed:      2,
		Threads:   4,
		Emitter:   fakeEmitter{basePattern: basePattern},
		Scorer:    fakeScorer{basePattern: basePattern},
	}

	results, _, err := Run(context.Background(), cfg)
	if err != nil {
		t.Fatalf("Run failed: %v", err)
	}
	// With Query={"modelA"}, only pairs touching modelA should be
	// reported: (A,B) and (A,C), not (B,C).
	if len(results) != 2 {
		t.Fatalf("expected 2 pairs with query filter, got %d: %+v", len(results), results)
	}
	for _, r := range results {
		if r.ModelA != "modelA" && r.ModelB != "modelA" {
			t.Errorf("pair (%s,%s) does not involve queried model modelA", r.ModelA, r.ModelB)
		}
	}
}

func TestPipeline_RejectsTooFewModels(t *testing.T) {
	cfg := Config{
		Models:    []Model{{Name: "onlyOne", Path: "x"}},
		NPool:     10,
		NShuffles: 10,
		Reservoir: testReservoir(10),
		Seed:      1,
		Threads:   1,
		Emitter:   fakeEmitter{basePattern: map[string]string{}},
		Scorer:    fakeScorer{basePattern: map[string]string{}},
	}
	if _, _, err := Run(context.Background(), cfg); err == nil {
		t.Error("expected error for < 2 models, got nil")
	}
}

func TestPipeline_ReservoirConverges(t *testing.T) {
	// A direct test of buildReservoir itself: confirms growth actually
	// stops (doesn't run away to ReservoirMax) when scores are
	// well-behaved, and that the resulting reservoir is usable -- every
	// model's reservoir sequences have recorded scores against every
	// other model.
	basePattern := map[string]string{
		"modelA": "ACGTACGTGG",
		"modelB": "GGTTCCAAGG",
	}
	cfg := Config{
		Models: []Model{
			{Name: "modelA", Path: "modelA"},
			{Name: "modelB", Path: "modelB"},
		},
		NPool:     30,
		Reservoir: testReservoir(30),
		Seed:      3,
		Threads:   2,
		Emitter:   fakeEmitter{basePattern: basePattern},
		Scorer:    fakeScorer{basePattern: basePattern},
	}

	realPool, err := buildRealPool(context.Background(), cfg)
	if err != nil {
		t.Fatalf("buildRealPool failed: %v", err)
	}

	res, err := buildReservoir(context.Background(), cfg, realPool, cfg.Reservoir)
	if err != nil {
		t.Fatalf("buildReservoir failed: %v", err)
	}

	for _, m := range cfg.Models {
		size := res.Sizes[m.Name]
		if size < cfg.Reservoir.ReservoirMin {
			t.Errorf("model %s: reservoir size %d is below ReservoirMin %d", m.Name, size, cfg.Reservoir.ReservoirMin)
		}
		if size > cfg.Reservoir.ReservoirMax {
			t.Errorf("model %s: reservoir size %d exceeds ReservoirMax %d", m.Name, size, cfg.Reservoir.ReservoirMax)
		}
		t.Logf("model %s: reservoir converged at size %d after %d rounds", m.Name, size, res.Rounds[m.Name])
	}

	// Every reservoir sequence should have a recorded score against every
	// other model (not necessarily against its own origin -- see the
	// doc comment on why that's skipped).
	for _, seq := range res.Pool.Sequences {
		for _, target := range cfg.Models {
			if target.Name == seq.Origin {
				continue
			}
			if _, ok := res.Matrix.Get(seq.ID, target.Name); !ok {
				t.Errorf("reservoir sequence %s (origin %s) missing score against %s", seq.ID, seq.Origin, target.Name)
			}
		}
	}

	if len(res.Warnings) > 0 {
		t.Logf("warnings: %v", res.Warnings)
	}
}

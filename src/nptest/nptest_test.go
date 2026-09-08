package nptest

import (
	"math"
	"testing"
)

func approxEqual(a, b, tol float64) bool { return math.Abs(a-b) <= tol }

func TestTwoSampleKS_IdenticalSamples(t *testing.T) {
	a := []float64{1, 2, 3, 4, 5}
	b := []float64{1, 2, 3, 4, 5}
	res := TwoSampleKS(a, b)
	if !approxEqual(res.D, 0, 1e-9) {
		t.Errorf("expected D=0 for identical samples, got %f", res.D)
	}
	if !approxEqual(res.PValue, 1.0, 1e-9) {
		t.Errorf("expected p=1 for identical samples, got %f", res.PValue)
	}
}

func TestTwoSampleKS_CompletelySeparated(t *testing.T) {
	a := []float64{1, 2, 3, 4, 5}
	b := []float64{100, 101, 102, 103, 104}
	res := TwoSampleKS(a, b)
	if !approxEqual(res.D, 1.0, 1e-9) {
		t.Errorf("expected D=1 for completely separated samples, got %f", res.D)
	}
	if res.PValue > 0.05 {
		t.Errorf("expected small p-value for completely separated samples, got %f", res.PValue)
	}
}

func TestTwoSampleKS_KnownHandComputedCase(t *testing.T) {
	// a = [1,2,3,4], b = [2,3,4,5]
	// F_a: at x=1 -> 1/4, x=2 -> 2/4, x=3 -> 3/4, x=4 -> 4/4, x=5 -> 4/4
	// F_b: at x=1 -> 0/4, x=2 -> 1/4, x=3 -> 2/4, x=4 -> 3/4, x=5 -> 4/4
	// max |diff| = 0.25
	a := []float64{1, 2, 3, 4}
	b := []float64{2, 3, 4, 5}
	res := TwoSampleKS(a, b)
	if !approxEqual(res.D, 0.25, 1e-9) {
		t.Errorf("expected D=0.25, got %f", res.D)
	}
}

func TestTwoSampleKS_PValueDecreasesWithSampleSize(t *testing.T) {
	aSmall := []float64{1, 2, 3}
	bSmall := []float64{2, 3, 4}
	aBig := make([]float64, 0, 300)
	bBig := make([]float64, 0, 300)
	for i := 0; i < 100; i++ {
		aBig = append(aBig, 1, 2, 3)
		bBig = append(bBig, 2, 3, 4)
	}
	small := TwoSampleKS(aSmall, bSmall)
	big := TwoSampleKS(aBig, bBig)
	if !approxEqual(small.D, big.D, 1e-9) {
		t.Errorf("expected same D regardless of replication, got small=%f big=%f", small.D, big.D)
	}
	if big.PValue >= small.PValue {
		t.Errorf("expected larger sample to give a smaller p-value for the same D, got small p=%f big p=%f", small.PValue, big.PValue)
	}
}

func TestMannWhitneyU_CompleteSeparationAHigher(t *testing.T) {
	a := []float64{10, 11, 12}
	b := []float64{1, 2, 3}
	res := MannWhitneyU(a, b)
	wantU := float64(len(a) * len(b))
	if !approxEqual(res.U, wantU, 1e-9) {
		t.Errorf("expected U=%f (complete separation), got %f", wantU, res.U)
	}
	if !approxEqual(res.AUC, 1.0, 1e-9) {
		t.Errorf("expected AUC=1.0, got %f", res.AUC)
	}
}

func TestMannWhitneyU_CompleteSeparationBHigher(t *testing.T) {
	a := []float64{1, 2, 3}
	b := []float64{10, 11, 12}
	res := MannWhitneyU(a, b)
	if !approxEqual(res.U, 0, 1e-9) {
		t.Errorf("expected U=0 when a is entirely lower, got %f", res.U)
	}
	if !approxEqual(res.AUC, 0.0, 1e-9) {
		t.Errorf("expected AUC=0.0, got %f", res.AUC)
	}
}

func TestMannWhitneyU_IdenticalInterleavedGivesAUCNearHalf(t *testing.T) {
	a := []float64{1, 3, 5, 7, 9}
	b := []float64{2, 4, 6, 8, 10}
	res := MannWhitneyU(a, b)
	if !approxEqual(res.AUC, 0.5, 0.15) {
		t.Errorf("expected AUC near 0.5 for interleaved samples, got %f", res.AUC)
	}
	if res.PValue < 0.5 {
		t.Errorf("expected a large (non-significant) p-value for evenly interleaved samples, got %f", res.PValue)
	}
}

func TestMannWhitneyU_KnownHandComputedCase(t *testing.T) {
	// a = [1,2,3], b = [4,5,6]: every a < every b.
	// U1 = rankSumA - n1(n1+1)/2 = (1+2+3) - 3*4/2 = 6 - 6 = 0
	a := []float64{1, 2, 3}
	b := []float64{4, 5, 6}
	res := MannWhitneyU(a, b)
	if !approxEqual(res.U, 0, 1e-9) {
		t.Errorf("expected U=0, got %f", res.U)
	}
}

func TestMannWhitneyU_TiesHandledSymmetrically(t *testing.T) {
	a := []float64{5, 5, 5}
	b := []float64{5, 5, 5}
	res := MannWhitneyU(a, b)
	if !approxEqual(res.AUC, 0.5, 1e-9) {
		t.Errorf("expected AUC=0.5 for fully tied samples, got %f", res.AUC)
	}
}

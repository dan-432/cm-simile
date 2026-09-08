// Package nptest implements the nonparametric two-sample tests used by
// cm-simile's KS+MW alternative run type (see src/ksmwcompare.go):
// Kolmogorov-Smirnov (any distributional difference) and Mann-Whitney U
// / Wilcoxon rank-sum (specifically powered against location shifts).
// Both operate on raw score samples directly -- unlike
// bitscoresumdifference's bootstrap-resampled-mean null, there is no
// resampling here, so these are only meaningful when applied to two
// genuine samples of individual observations (e.g. a real pool's scores
// vs. a shuffled-null pool's scores), not to a distribution of already-
// averaged statistics.
package nptest

import (
	"math"
	"sort"
)

// KSResult holds the result of a two-sample Kolmogorov-Smirnov test.
type KSResult struct {
	D      float64 // max absolute difference between empirical CDFs, in [0,1]
	PValue float64 // asymptotic two-sided p-value
	N1, N2 int
}

// TwoSampleKS performs the two-sample Kolmogorov-Smirnov test comparing
// the empirical distributions of a and b, using the standard asymptotic
// (Kolmogorov) approximation for the p-value -- the same approximation
// most standard implementations (e.g. SciPy, R's ks.test for
// moderate-to-large samples) use.
//
// KS tests for ANY distributional difference (location, spread, shape),
// not just a shift in central tendency. This makes it a useful
// complement to MannWhitneyU (which is specifically powered against
// shifts): a pair whose real pool is a MIXTURE -- some sequences
// retaining strong family-specific signal, others looking noise-like,
// e.g. from partial/degraded homology -- can produce a bimodal real
// distribution against a unimodal null that a location-shift test
// underpowers but KS's max-CDF-gap statistic can still catch.
//
// D itself, not just the p-value, is worth reporting directly: at
// reservoir/large-sample scale, KS p-values become extremely sensitive
// (effective sample size n1*n2/(n1+n2) can be large even for modest
// per-sample sizes), so a "significant" p-value can reflect a tiny,
// unimportant distributional wobble. D is bounded in [0,1] and
// interpretable on its own as "the largest gap between the two ECDFs",
// independent of sample size.
func TwoSampleKS(a, b []float64) KSResult {
	n1, n2 := len(a), len(b)

	sortedA := append([]float64(nil), a...)
	sortedB := append([]float64(nil), b...)
	sort.Float64s(sortedA)
	sort.Float64s(sortedB)

	combined := make([]float64, 0, n1+n2)
	combined = append(combined, sortedA...)
	combined = append(combined, sortedB...)
	sort.Float64s(combined)

	var d float64
	for _, x := range combined {
		f1 := float64(countLE(sortedA, x)) / float64(n1)
		f2 := float64(countLE(sortedB, x)) / float64(n2)
		if diff := math.Abs(f1 - f2); diff > d {
			d = diff
		}
	}

	ne := float64(n1*n2) / float64(n1+n2)
	p := ksAsymptoticPValue(math.Sqrt(ne) * d)

	return KSResult{D: d, PValue: p, N1: n1, N2: n2}
}

// countLE returns the number of elements in sorted (already ascending)
// that are <= x.
func countLE(sorted []float64, x float64) int {
	return sort.Search(len(sorted), func(i int) bool { return sorted[i] > x })
}

// ksAsymptoticPValue evaluates the Kolmogorov distribution's survival
// function Q(t) = 2 * sum_{k=1}^inf (-1)^(k-1) exp(-2 k^2 t^2), truncated
// once terms become numerically negligible (they decay very fast --
// this converges in well under 100 terms for any t of practical
// interest).
func ksAsymptoticPValue(t float64) float64 {
	if t < 1e-10 {
		return 1.0
	}
	var sum float64
	sign := 1.0
	for k := 1; k <= 100; k++ {
		term := sign * math.Exp(-2*float64(k*k)*t*t)
		sum += term
		sign = -sign
		if math.Abs(term) < 1e-12 {
			break
		}
	}
	p := 2 * sum
	if p < 0 {
		p = 0
	}
	if p > 1 {
		p = 1
	}
	return p
}

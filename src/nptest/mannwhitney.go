package nptest

import (
	"math"
	"sort"
)

// MannWhitneyResult holds the result of a two-sample Mann-Whitney U test
// (equivalently, the Wilcoxon rank-sum test).
type MannWhitneyResult struct {
	U      float64 // Mann-Whitney U statistic for sample a
	AUC    float64 // U / (n1*n2): P(random a-value > random b-value), ties counted as a half-win
	Z      float64 // tie-corrected normal-approximation z-score
	PValue float64 // two-sided p-value from the normal approximation
	N1, N2 int
}

// MannWhitneyU performs the two-sample Mann-Whitney U test, specifically
// powered against location shifts (does one sample tend to have higher
// values than the other) -- see TwoSampleKS for a test sensitive to any
// distributional difference, including shape differences a pure shift
// test can miss.
//
// Implemented directly (rank combined samples, tie-corrected normal
// approximation) rather than via an external stats library, matching
// the same approach already used for Mann-Whitney/AUC in
// scripts/run_clan_benchmark.py -- this is the Go port of that same
// well-tested logic, so the two should agree on identical inputs.
func MannWhitneyU(a, b []float64) MannWhitneyResult {
	n1, n2 := len(a), len(b)

	type tagged struct {
		val   float64
		group int // 0 = a, 1 = b
	}
	combined := make([]tagged, 0, n1+n2)
	for _, v := range a {
		combined = append(combined, tagged{v, 0})
	}
	for _, v := range b {
		combined = append(combined, tagged{v, 1})
	}
	sort.Slice(combined, func(i, j int) bool { return combined[i].val < combined[j].val })

	ranks := make([]float64, len(combined))
	var tieTerm float64
	i := 0
	for i < len(combined) {
		j := i
		for j < len(combined) && combined[j].val == combined[i].val {
			j++
		}
		avgRank := float64(i+1+j) / 2.0 // average of ranks i+1..j, 1-indexed
		for k := i; k < j; k++ {
			ranks[k] = avgRank
		}
		t := float64(j - i)
		tieTerm += t*t*t - t
		i = j
	}

	var rankSumA float64
	for idx, c := range combined {
		if c.group == 0 {
			rankSumA += ranks[idx]
		}
	}
	u1 := rankSumA - float64(n1)*float64(n1+1)/2.0
	var auc float64
	if n1 > 0 && n2 > 0 {
		auc = u1 / (float64(n1) * float64(n2))
	}

	n := float64(n1 + n2)
	meanU := float64(n1) * float64(n2) / 2.0
	var varU float64
	if n > 1 {
		varU = (float64(n1) * float64(n2) / 12.0) * ((n + 1) - tieTerm/(n*(n-1)))
	}

	z, p := math.NaN(), math.NaN()
	if varU > 0 {
		z = (u1 - meanU) / math.Sqrt(varU)
		p = 2.0 * (1.0 - normCDF(math.Abs(z)))
	}

	return MannWhitneyResult{U: u1, AUC: auc, Z: z, PValue: p, N1: n1, N2: n2}
}

func normCDF(x float64) float64 {
	return 0.5 * (1.0 + math.Erf(x/math.Sqrt2))
}

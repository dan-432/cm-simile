// Package stats computes the directional overlap statistics and bootstrap
// p-values that are cm-simile's actual output, given a real score matrix
// and a set of null (dinucleotide-shuffled) score matrices already
// populated by the pipeline package.
package stats

import (
	"fmt"

	"cmsimile/src/matrix"
	"cmsimile/src/seqpool"
)

// PairStat is the full result of comparing two models, A and B.
//
// The two directional statistics answer different questions and are kept
// separate rather than collapsed into one symmetric number:
//   - StatAintoB = mean score of A's real sequences against B, minus mean
//     score of A's real sequences against A itself. Near zero means B
//     recognises A's typical output about as well as A does; very
//     negative means B does not recognise A's sequences.
//   - StatBintoA is the same computation in the other direction.
//
// A pair where both directions are near zero AND both p-values are small
// (i.e. the observed cross-recognition is unlikely under the
// composition-matched null) indicates the two models are substantially
// redundant. A pair significant in only one direction indicates
// asymmetric containment: one model's typical sequences are recognised by
// the other, but not vice versa (a signature of a broader model
// containing a narrower one, e.g. a family within a more permissive
// superfamily-level model).
type PairStat struct {
	ModelA, ModelB string

	// MeanScoreAonA is A's real pool scored against A (i.e. against
	// itself); MeanScoreAonB is A's real pool scored against B; and so on.
	MeanScoreAonA, MeanScoreAonB float64
	MeanScoreBonB, MeanScoreBonA float64

	StatAintoB float64 // mean(A on B) - mean(A on A)
	StatBintoA float64 // mean(B on A) - mean(B on B)

	// PValueAintoB is the fraction of null replicates whose shuffled-A
	// scored-on-B statistic is >= the observed StatAintoB (one-sided:
	// "cross-recognition at least this strong"). Likewise PValueBintoA.
	PValueAintoB float64
	PValueBintoA float64

	NShuffles int
}

// ComputePairStat computes the full PairStat for models A and B, given:
//   - real: the score matrix built from the real (unshuffled) pool,
//     containing scores for every model against every real sequence.
//   - nullAonB: for each shuffle replicate, the mean score of shuffled-A
//     sequences against model B (pre-aggregated by the pipeline, since the
//     per-sequence shuffled scores themselves aren't needed downstream).
//   - nullBonA: likewise, mean score of shuffled-B sequences against
//     model A, per replicate.
//
// pool must contain the real (unshuffled) sequences tagged with origins A
// and B, so their IDs can be looked up in `real`.
func ComputePairStat(pool *seqpool.Pool, real *matrix.ScoreMatrix, modelA, modelB string, nullAonB, nullBonA []float64) (PairStat, error) {
	idsA := ids(pool.ByOrigin(modelA))
	idsB := ids(pool.ByOrigin(modelB))
	if len(idsA) == 0 {
		return PairStat{}, fmt.Errorf("stats: no real pool sequences found for model %q", modelA)
	}
	if len(idsB) == 0 {
		return PairStat{}, fmt.Errorf("stats: no real pool sequences found for model %q", modelB)
	}

	meanAonA, err := real.MeanFor(idsA, modelA)
	if err != nil {
		return PairStat{}, err
	}
	meanAonB, err := real.MeanFor(idsA, modelB)
	if err != nil {
		return PairStat{}, err
	}
	meanBonB, err := real.MeanFor(idsB, modelB)
	if err != nil {
		return PairStat{}, err
	}
	meanBonA, err := real.MeanFor(idsB, modelA)
	if err != nil {
		return PairStat{}, err
	}

	statAintoB := meanAonB - meanAonA
	statBintoA := meanBonA - meanBonB

	if len(nullAonB) != len(nullBonA) {
		return PairStat{}, fmt.Errorf("stats: mismatched null replicate counts (%d vs %d)", len(nullAonB), len(nullBonA))
	}

	pAintoB := empiricalPValue(statAintoB, nullAonB, meanAonA)
	pBintoA := empiricalPValue(statBintoA, nullBonA, meanBonB)

	return PairStat{
		ModelA: modelA, ModelB: modelB,
		MeanScoreAonA: meanAonA, MeanScoreAonB: meanAonB,
		MeanScoreBonB: meanBonB, MeanScoreBonA: meanBonA,
		StatAintoB: statAintoB, StatBintoA: statBintoA,
		PValueAintoB: pAintoB, PValueBintoA: pBintoA,
		NShuffles: len(nullAonB),
	}, nil
}

// empiricalPValue computes the one-sided empirical p-value: the fraction
// of null replicate statistics (each: mean shuffled-cross-score minus the
// real self-score baseline) that are >= the observed statistic. Using the
// same real self-score baseline for both observed and null statistics
// keeps them on the same scale, since only the "cross" side of the
// subtraction differs between the real and shuffled cases.
//
// nullCrossScores are per-replicate mean scores of the shuffled pool
// against the *other* model (e.g. shuffled-A scored against B); selfBaseline
// is the real self-score mean (e.g. real-A scored against A) used to form
// each null replicate's statistic the same way the observed statistic was
// formed.
func empiricalPValue(observedStat float64, nullCrossScores []float64, selfBaseline float64) float64 {
	if len(nullCrossScores) == 0 {
		return 1.0 // no null replicates: cannot reject, maximally conservative
	}
	countAtLeastAsExtreme := 0
	for _, crossScore := range nullCrossScores {
		nullStat := crossScore - selfBaseline
		if nullStat >= observedStat {
			countAtLeastAsExtreme++
		}
	}
	// +1/+1 continuity correction: a p-value of exactly 0 from finite
	// resampling overstates precision (see North, Curtis & Sham 2002);
	// this bounds the smallest reportable p-value at 1/(n+1) rather than
	// claiming impossible certainty from a finite bootstrap.
	return float64(countAtLeastAsExtreme+1) / float64(len(nullCrossScores)+1)
}

func ids(seqs []seqpool.Sequence) []string {
	out := make([]string, len(seqs))
	for i, s := range seqs {
		out[i] = s.ID
	}
	return out
}

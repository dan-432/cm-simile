// Package matrix holds the (sequence x model) score matrix that pairwise
// statistics are computed from. Building this once per pool -- rather than
// re-scoring per pair -- is what keeps cm-simile's subprocess call count
// linear in the number of models rather than quadratic in the number of
// pairs being compared.
package matrix

import (
	"sort"
	"fmt"
)

// ScoreMatrix holds bit scores for a pool of sequences scored against a
// set of covariance models: Scores[seqID][modelName] = bit score. It is
// deliberately just a map of maps rather than a dense 2D array, since in
// the bipartite (query-vs-reference) use case most (sequence, model) pairs
// are never scored at all -- a dense matrix would be mostly empty.
type ScoreMatrix struct {
	scores map[string]map[string]float64
}

// New returns an empty ScoreMatrix.
func New() *ScoreMatrix {
	return &ScoreMatrix{scores: make(map[string]map[string]float64)}
}

// Set records the bit score of sequence seqID against model.
func (m *ScoreMatrix) Set(seqID, model string, score float64) {
	if m.scores[seqID] == nil {
		m.scores[seqID] = make(map[string]float64)
	}
	m.scores[seqID][model] = score
}

// SetAll records every entry of scores (seqID -> bit score) against model
// in one call -- the usual way a ScoreMatrix gets populated, since a
// single cmscore invocation returns scores for a whole pool at once.
func (m *ScoreMatrix) SetAll(model string, scores map[string]float64) {
	for seqID, sc := range scores {
		m.Set(seqID, model, sc)
	}
}

// Get returns the bit score of seqID against model, and whether it was
// present.
func (m *ScoreMatrix) Get(seqID, model string) (float64, bool) {
	row, ok := m.scores[seqID]
	if !ok {
		return 0, false
	}
	sc, ok := row[model]
	return sc, ok
}

// MustGet is like Get but returns an error instead of a boolean, for
// call sites where a missing entry indicates a real bug (e.g. incomplete
// scoring) rather than an expected sparse gap.
func (m *ScoreMatrix) MustGet(seqID, model string) (float64, error) {
	sc, ok := m.Get(seqID, model)
	if !ok {
		return 0, fmt.Errorf("matrix: no score recorded for sequence %q against model %q", seqID, model)
	}
	return sc, nil
}

// Row returns a copy of the score map for one sequence -- model name ->
// bit score, for every model that sequence has been scored against.
// Returns nil if the sequence has no recorded scores. A copy (not the
// internal map) is returned so callers can't accidentally mutate matrix
// state through it. Used for dumping full score matrices (see
// pipeline.DumpMatrices) without needing a caller to already know which
// models a given sequence was scored against.
func (m *ScoreMatrix) Row(seqID string) map[string]float64 {
	row, ok := m.scores[seqID]
	if !ok {
		return nil
	}
	out := make(map[string]float64, len(row))
	for k, v := range row {
		out[k] = v
	}
	return out
}

// MeanFor computes the mean bit score, against model, of the sequences in
// seqIDs. Returns an error if any sequence lacks a score against model
// (a scoring gap should surface loudly, not silently skew the mean over
// fewer sequences than intended).
func (m *ScoreMatrix) MeanFor(seqIDs []string, model string) (float64, error) {
	if len(seqIDs) == 0 {
		return 0, fmt.Errorf("matrix: MeanFor called with no sequence IDs")
	}
	var sum float64
	for _, id := range seqIDs {
		sc, err := m.MustGet(id, model)
		if err != nil {
			return 0, err
		}
		sum += sc
	}
	return sum / float64(len(seqIDs)), nil
}

// MedianFor computes the median bit score, against model, of the sequences in
// seqIDs. Returns an error if any sequence lacks a score against model
// (a scoring gap should surface loudly, not silently skew the median over
// fewer sequences than intended).
func (m *ScoreMatrix) MedianFor(seqIDs []string, model string) (float64, error) {
	if len(seqIDs) == 0 {
		return 0, fmt.Errorf("matrix: MedianFor called with no sequence IDs")
	}
	var scores []float64
	for _, id := range seqIDs {
		sc, err := m.MustGet(id, model)
		if err != nil {
			return 0, err
		}
		scores = append(scores, sc)
	}
	sort.Float64s(scores)
	n := len(scores)
	if n%2 == 1 {
		return scores[n/2], nil
	}
	return (scores[n/2-1] + scores[n/2]) / 2, nil
}


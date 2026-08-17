// Package output writes cm-simile results in plain TSV, suitable for
// downstream analysis in R/pandas or direct inspection.
package output

import (
	"fmt"
	"io"

	"cmsimile/src/stats"
)

// WriteTSV writes one row per PairStat, tab-separated, with a header row.
func WriteTSV(w io.Writer, results []stats.PairStat) error {
	header := []string{
		"model_a", "model_b",
		"median_score_a_on_a", "median_score_a_on_b",
		"median_score_b_on_b", "median_score_b_on_a",
		"stat_a_into_b", "stat_b_into_a",
		"pvalue_a_into_b", "pvalue_b_into_a",
		"n_shuffles",
	}
	if _, err := fmt.Fprintln(w, joinTab(header)); err != nil {
		return err
	}
	for _, r := range results {
		row := []string{
			r.ModelA, r.ModelB,
			ftoa(r.MedianScoreAonA), ftoa(r.MedianScoreAonB),
			ftoa(r.MedianScoreBonB), ftoa(r.MedianScoreBonA),
			ftoa(r.StatAintoB), ftoa(r.StatBintoA),
			ftoa(r.PValueAintoB), ftoa(r.PValueBintoA),
			itoa(r.NShuffles),
		}
		if _, err := fmt.Fprintln(w, joinTab(row)); err != nil {
			return err
		}
	}
	return nil
}

func joinTab(fields []string) string {
	out := fields[0]
	for _, f := range fields[1:] {
		out += "\t" + f
	}
	return out
}

func ftoa(f float64) string { return fmt.Sprintf("%.4f", f) }
func itoa(i int) string     { return fmt.Sprintf("%d", i) }

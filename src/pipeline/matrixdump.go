package pipeline

import (
	"fmt"
	"os"
	"path/filepath"
	"sort"

	"cmsimile/src/matrix"
	"cmsimile/src/seqpool"
)

// DumpMatrices writes the real pool's and reservoir's full (sequence x
// model) score matrices to dir as long-format TSVs: one row per
// (sequence, model) pair actually scored -- sequence_id, origin, model,
// score. Long format rather than wide (one column per model) since a run
// can involve many models; long format needs no ragged-column handling
// and loads directly into pandas (`pd.read_csv(..., sep="\t")`) or R
// (`readr::read_tsv`), with a pivot to wide left to the caller if they
// want it.
//
// reservoir_matrix.tsv in particular is what you want for plotting the
// null distribution's actual shape (e.g. a histogram of one model's
// reservoir scores against a specific target) rather than only the
// summary statistics in the main results TSV.
func DumpMatrices(dir string, realPool *seqpool.Pool, realMatrix *matrix.ScoreMatrix, reservoirPool *seqpool.Pool, reservoirMatrix *matrix.ScoreMatrix) error {
	if err := os.MkdirAll(dir, 0o755); err != nil {
		return fmt.Errorf("creating matrix dump directory: %w", err)
	}
	if err := dumpOne(filepath.Join(dir, "real_matrix.tsv"), realPool, realMatrix); err != nil {
		return fmt.Errorf("dumping real matrix: %w", err)
	}
	if err := dumpOne(filepath.Join(dir, "reservoir_matrix.tsv"), reservoirPool, reservoirMatrix); err != nil {
		return fmt.Errorf("dumping reservoir matrix: %w", err)
	}
	return nil
}

func dumpOne(path string, pool *seqpool.Pool, m *matrix.ScoreMatrix) error {
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()

	if _, err := fmt.Fprintln(f, "sequence_id\torigin\tmodel\tscore"); err != nil {
		return err
	}

	for _, seq := range pool.Sequences {
		row := m.Row(seq.ID)
		if row == nil {
			continue // sequence was pooled but never scored (shouldn't normally happen); skip rather than error, since a dump is diagnostic, not load-bearing
		}
		// Sort model names for deterministic output -- map iteration
		// order would otherwise make every dump differ byte-for-byte
		// between runs even with identical data, which is a minor but
		// real annoyance for diffing/version-controlling dumps.
		models := make([]string, 0, len(row))
		for model := range row {
			models = append(models, model)
		}
		sort.Strings(models)

		for _, model := range models {
			if _, err := fmt.Fprintf(f, "%s\t%s\t%s\t%.6f\n", seq.ID, seq.Origin, model, row[model]); err != nil {
				return err
			}
		}
	}
	return nil
}

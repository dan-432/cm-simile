package infernal

import (
	"bufio"
	"context"
	"fmt"
	"os"
	"strconv"
	"strings"

	"cmsimile/src/seqpool"
)

// Scorer scores a pool of sequences against a covariance model using
// Infernal's cmsearch, run with all acceleration filters disabled and an
// inclusion threshold low enough that every sequence in the pool is
// guaranteed to appear in the output -- turning what is normally a
// "find significant hits" tool into a "score every sequence, however
// weak" tool, which is what cm-simile's statistics need.
//
// This replaces an earlier cmscore-based implementation. cmscore's plain
// text tabular output format isn't documented with a stable column
// layout and has shifted across Infernal releases, whereas cmsearch's
// --tblout format is a fixed, documented 18-column layout that doesn't
// change release to release -- see http://eddylab.org/infernal/Userguide.pdf,
// "tabular output" -- making it the more reliable thing to parse
// programmatically.
type Scorer struct {
	// BinPath is the path to the cmsearch binary. If empty, "cmsearch" is
	// resolved from PATH at call time.
	BinPath string

	// CPU is passed as cmsearch's --cpu, controlling cmsearch's own
	// internal multithreading for a single invocation. This is separate
	// from (and multiplies with) the pipeline's own --threads, which
	// controls how many concurrent cmsearch subprocesses are dispatched
	// at once -- see pipeline.Config.Threads. Defaults to 1 if unset.
	CPU int

	// MinBitScore is passed as cmsearch's -T, a raw bit-score inclusion
	// threshold. cm-simile needs a real score for every sequence in the
	// pool -- including dinucleotide-shuffled null sequences that are
	// *expected* to score badly -- so this should be set low enough that
	// no sequence's best hit is ever excluded from --tblout by score
	// alone. Using -T (bit score) rather than -E (E-value) avoids
	// E-value's dependence on database/pool size shifting the effective
	// cutoff between the real-pool and null-pool scoring passes. Defaults
	// to -1000.0 if zero (Go's zero value), which is low enough for any
	// realistically sized CM.
	MinBitScore float64

	// DisableFilters, if true (the default), passes --max to cmsearch,
	// disabling all HMM-based acceleration filtering so every sequence is
	// scored by the full CM. Without this, cmsearch's default filter
	// pipeline can discard a weak-scoring sequence before it ever reaches
	// the CM scoring stage -- which would silently drop exactly the
	// low-scoring shuffled-null sequences cm-simile's p-values depend on
	// having a real score for. This makes cmsearch significantly slower
	// (full Inside/CYK on every sequence, no shortcuts), which is an
	// accepted cost for correctness at cm-simile's pool sizes (tens to
	// low hundreds of sequences per model, not a genome-scale database).
	DisableFilters bool

	// KeepTempFiles, if true, leaves the temporary FASTA and --tblout
	// files on disk instead of deleting them (useful for debugging a
	// parsing mismatch by hand against real cmsearch output).
	KeepTempFiles bool
}

func (s Scorer) bin() string { return resolveBin(s.BinPath, "cmsearch") }

func (s Scorer) cpu() int {
	if s.CPU <= 0 {
		return 1
	}
	return s.CPU
}

func (s Scorer) minBitScore() float64 {
	if s.MinBitScore == 0 {
		return -1000.0
	}
	return s.MinBitScore
}

// Score runs cmsearch against every sequence in seqs and returns a map
// from sequence ID to its best (maximum) bit score against cmPath.
// "Best" matters here because cmsearch reports one row per hit *region*,
// not per sequence -- a sequence can have multiple weak partial hits in
// different frames/positions, and the maximum is the representative score
// for "how well does this model recognise this sequence overall".
//
// Command run (see the DisableFilters/MinBitScore doc comments for why
// these specific flags):
//
//	cmsearch --cpu <CPU> -T <MinBitScore> --max --noali --tblout <tmp> <cmPath> <tmpFasta>
func (s Scorer) Score(ctx context.Context, cmPath string, seqs []seqpool.Sequence) (map[string]float64, error) {
	if len(seqs) == 0 {
		return nil, fmt.Errorf("infernal: Score called with empty sequence set")
	}

	tmpFasta, err := os.CreateTemp("", "cmsimile-pool-*.fasta")
	if err != nil {
		return nil, fmt.Errorf("cmsearch: creating temp FASTA: %w", err)
	}
	fastaPath := tmpFasta.Name()
	defer func() {
		if !s.KeepTempFiles {
			os.Remove(fastaPath)
		}
	}()

	pool := seqpool.Pool{Sequences: seqs}
	if err := pool.WriteFASTA(tmpFasta); err != nil {
		tmpFasta.Close()
		return nil, fmt.Errorf("cmsearch: writing temp FASTA: %w", err)
	}
	if err := tmpFasta.Close(); err != nil {
		return nil, fmt.Errorf("cmsearch: closing temp FASTA: %w", err)
	}

	tblFile, err := os.CreateTemp("", "cmsimile-tblout-*.tbl")
	if err != nil {
		return nil, fmt.Errorf("cmsearch: creating temp tblout file: %w", err)
	}
	tblPath := tblFile.Name()
	tblFile.Close()
	defer func() {
		if !s.KeepTempFiles {
			os.Remove(tblPath)
		}
	}()

	args := []string{
		"--cpu", strconv.Itoa(s.cpu()),
		"-T", strconv.FormatFloat(s.minBitScore(), 'f', -1, 64),
		"--noali",
		"--tblout", tblPath,
	}
	if s.DisableFilters {
		args = append(args, "--max")
	}
	args = append(args, cmPath, fastaPath)

	if _, err := runCapture(ctx, s.bin(), args...); err != nil {
		return nil, fmt.Errorf("cmsearch against %s: %w", cmPath, err)
	}

	scores, err := parseTblout(tblPath)
	if err != nil {
		return nil, fmt.Errorf("cmsearch against %s: parsing --tblout: %w", cmPath, err)
	}

	for _, sq := range seqs {
		if _, ok := scores[sq.ID]; !ok {
			return nil, fmt.Errorf(
				"cmsearch against %s: no hit reported for sequence %q -- if DisableFilters is false this may mean the acceleration filter dropped a weak sequence; if MinBitScore is set high this may mean a real hit scored below it",
				cmPath, sq.ID)
		}
	}
	return scores, nil
}

// tbloutScoreColumn / tbloutNameColumn: 0-indexed positions of the target
// (sequence) name and bit score columns in cmsearch's --tblout format.
// This layout is fixed and documented (18 whitespace-delimited columns,
// the last -- description -- free-text and irrelevant here), unlike
// cmscore's plain-text table, which is why cm-simile scores via cmsearch
// rather than cmscore. See the Infernal user guide's "tabular output"
// section for the full column reference.
const (
	tbloutNameColumn  = 0 // target name -- the sequence we queried
	tbloutScoreColumn = 14
)

// parseTblout parses cmsearch's --tblout format, returning the maximum
// bit score seen per target (sequence) name across all reported hit rows.
func parseTblout(path string) (map[string]float64, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	scores := make(map[string]float64)
	scanner := bufio.NewScanner(f)
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		fields := strings.Fields(line)
		if len(fields) <= tbloutScoreColumn {
			continue // malformed/short line, not a real data row
		}
		name := fields[tbloutNameColumn]
		score, err := strconv.ParseFloat(fields[tbloutScoreColumn], 64)
		if err != nil {
			continue
		}
		// A sequence can have multiple hit regions reported as separate
		// rows; keep the best (maximum) as the representative score.
		if existing, ok := scores[name]; !ok || score > existing {
			scores[name] = score
		}
	}
	if err := scanner.Err(); err != nil {
		return nil, err
	}
	return scores, nil
}

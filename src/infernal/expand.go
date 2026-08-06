package infernal

import (
	"bufio"
	"context"
	"fmt"
	"os"
	"path/filepath"
	"strings"
)

// ExpandedModel names one covariance model that has been materialized as
// its own single-model file on disk, ready to pass to cmemit/cmsearch
// (both of which expect a single-model CM file, not an index into a
// multi-model one).
type ExpandedModel struct {
	Name string // the model's name as recorded in the CM file (via cmstat)
	Path string // path to the single-model .cm file on disk
}

// Expander materializes individual models out of a (possibly multi-model)
// CM file, via cmstat (to list model names) and cmfetch (to extract each
// one). This is what lets cm-simile take `Rfam.cm` directly as input
// rather than requiring it pre-split into one file per family.
type Expander struct {
	CmstatBin  string
	CmfetchBin string
}

func (e Expander) cmstat() string  { return resolveBin(e.CmstatBin, "cmstat") }
func (e Expander) cmfetch() string { return resolveBin(e.CmfetchBin, "cmfetch") }

// ExpandCMFile expands cmPath into one ExpandedModel per named model it
// contains, materializing each into its own file under workDir.
//
//   - If cmPath contains a single model, one ExpandedModel is returned,
//     backed by a copy of cmPath itself (no cmfetch needed).
//   - If cmPath contains many models (e.g. all of Rfam), each is
//     extracted via `cmfetch <cmPath> <name>` into workDir. For files
//     with many models, an SSI index is built once up front
//     (`cmfetch --index`) if one doesn't already exist alongside cmPath,
//     so each individual fetch is a fast indexed lookup rather than a
//     linear scan of the whole file -- the difference between this being
//     fast and this being O(N^2) at Rfam scale (~4,100+ families).
func (e Expander) ExpandCMFile(ctx context.Context, cmPath, workDir string) ([]ExpandedModel, error) {
	names, err := e.listModelNames(ctx, cmPath)
	if err != nil {
		return nil, err
	}

	if len(names) <= 1 {
		// Single (or unnamed/unparseable-as-multi) model: just use the
		// file as-is under its base filename, no cmfetch round-trip
		// needed.
		name := names[0]
		if name == "" {
			name = strings.TrimSuffix(filepath.Base(cmPath), filepath.Ext(cmPath))
		}
		dst := filepath.Join(workDir, sanitizeFilename(name)+".cm")
		if err := copyFile(cmPath, dst); err != nil {
			return nil, fmt.Errorf("infernal: copying single-model CM %s: %w", cmPath, err)
		}
		return []ExpandedModel{{Name: name, Path: dst}}, nil
	}

	if err := e.ensureSSIIndex(ctx, cmPath); err != nil {
		return nil, fmt.Errorf("infernal: building SSI index for %s: %w", cmPath, err)
	}

	out := make([]ExpandedModel, 0, len(names))
	for _, name := range names {
		dst := filepath.Join(workDir, sanitizeFilename(name)+".cm")
		data, err := runCapture(ctx, e.cmfetch(), cmPath, name)
		if err != nil {
			return nil, fmt.Errorf("infernal: cmfetch %s from %s: %w", name, cmPath, err)
		}
		if err := os.WriteFile(dst, data, 0o644); err != nil {
			return nil, fmt.Errorf("infernal: writing extracted model %s: %w", name, err)
		}
		out = append(out, ExpandedModel{Name: name, Path: dst})
	}
	return out, nil
}

// ensureSSIIndex builds a cmfetch SSI index for cmPath (cmPath + ".ssi")
// if one doesn't already exist. Building it once up front turns each
// subsequent cmfetch call into an O(1) indexed lookup instead of an O(N)
// scan of the whole file -- essential for expanding something the size of
// Rfam.cm without quadratic blowup.
func (e Expander) ensureSSIIndex(ctx context.Context, cmPath string) error {
	if _, err := os.Stat(cmPath + ".ssi"); err == nil {
		return nil // already indexed
	}
	_, err := runCapture(ctx, e.cmfetch(), "--index", cmPath)
	return err
}

// listModelNames runs `cmstat <cmPath>` and extracts the model name from
// each data row of its summary table. cmstat's table format has an idx
// column followed by a name column; header/comment lines are skipped.
func (e Expander) listModelNames(ctx context.Context, cmPath string) ([]string, error) {
	out, err := runCapture(ctx, e.cmstat(), cmPath)
	if err != nil {
		return nil, fmt.Errorf("cmstat failed for %s: %w", cmPath, err)
	}

	seen := make(map[string]bool)
	var names []string
	scanner := bufio.NewScanner(strings.NewReader(string(out)))
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		fields := strings.Fields(line)
		if len(fields) < 2 {
			continue
		}
		name := fields[1]
		if name == "name" || name == "idx" { // stray header-ish tokens
			continue
		}
		if seen[name] {
			continue
		}
		seen[name] = true
		names = append(names, name)
	}
	if err := scanner.Err(); err != nil {
		return nil, err
	}
	if len(names) == 0 {
		// cmstat produced no parseable rows -- treat as a single
		// (possibly unnamed) model rather than failing outright.
		return []string{""}, nil
	}
	return names, nil
}

func sanitizeFilename(name string) string {
	replacer := strings.NewReplacer(" ", "_", "/", "_", ":", "_", "\\", "_", "\t", "_")
	return replacer.Replace(strings.TrimSpace(name))
}

func copyFile(src, dst string) error {
	data, err := os.ReadFile(src)
	if err != nil {
		return err
	}
	return os.WriteFile(dst, data, 0o644)
}

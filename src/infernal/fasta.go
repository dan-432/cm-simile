package infernal

import (
	"bufio"
	"fmt"
	"io"
	"strings"
)

// FastaRecord is a single parsed FASTA entry: the header line (without the
// leading '>') and the concatenated, line-unwrapped sequence.
type FastaRecord struct {
	Header string
	Seq    string
}

// ParseFASTA parses simple FASTA text (as produced by cmemit) into a slice
// of records. It tolerates wrapped sequence lines and blank lines, and does
// not validate alphabet -- callers that care about a specific alphabet
// (e.g. seqpool.NewSequence normalising U->T) should do so downstream.
func ParseFASTA(r io.Reader) ([]FastaRecord, error) {
	scanner := bufio.NewScanner(r)
	// cmemit output for long sampled sequences can exceed the default
	// bufio.Scanner token size on a single wrapped line in unusual cases;
	// give it plenty of headroom.
	buf := make([]byte, 0, 64*1024)
	scanner.Buffer(buf, 16*1024*1024)

	var records []FastaRecord
	var cur *FastaRecord
	var seqBuilder strings.Builder

	flush := func() {
		if cur != nil {
			cur.Seq = seqBuilder.String()
			records = append(records, *cur)
		}
	}

	for scanner.Scan() {
		line := strings.TrimRight(scanner.Text(), "\r\n")
		if line == "" {
			continue
		}
		if strings.HasPrefix(line, ">") {
			flush()
			header := strings.TrimSpace(strings.TrimPrefix(line, ">"))
			cur = &FastaRecord{Header: header}
			seqBuilder.Reset()
			continue
		}
		if cur == nil {
			return nil, fmt.Errorf("infernal: FASTA data before any header line: %q", line)
		}
		seqBuilder.WriteString(strings.TrimSpace(line))
	}
	flush()

	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("infernal: reading FASTA: %w", err)
	}
	return records, nil
}

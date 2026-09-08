package mafft

import (
	"context"
	"fmt"
	"os"
	"os/exec"
)

// runCapture runs bin with args, piping stdout to the specified output file. stderr is
// included in the returned error on failure so a failing Mafft
// invocation is diagnosable from the error
// message alone, without needing to re-run by hand.
func runCapture(ctx context.Context, bin string, output string, args ...string) error {
	file, err := os.Create(output)
	if err != nil {
	    return fmt.Errorf("failed to create output file: %w", err)
	}
	defer file.Close()

	// Prepare the command (e.g., "echo hello")
	cmd := exec.CommandContext(ctx, bin, args...)

	// Redirect stdout directly to the file
	cmd.Stdout = file

	// Run the command
	err = cmd.Run()
	if err != nil {
	    return fmt.Errorf("command failed: %w", err)
	}

	return nil
}


// func runCapture(ctx context.Context, bin string, output string, args ...string) error {
// 	cmd := exec.CommandContext(ctx, bin, args...)

// 	// Create a pipe to read stdout
// 	stdout, err := cmd.StdoutPipe()
// 	if err != nil {
// 		return fmt.Errorf("failed to create stdout pipe: %w", err)
// 	}

// 	// Start the command asynchronously
// 	if err := cmd.Start(); err != nil {
// 		return fmt.Errorf("failed to start command: %w", err)
// 	}

// 	// Stream and process the output line by line
// 	scanner := bufio.NewScanner(stdout)
// 	for scanner.Scan() {
// 		line := scanner.Text()
// 		// Write the line to the output file
		
// 	}

// 	// Wait for the command to finish after the pipe closes
// 	if err := cmd.Wait(); err != nil {
// 		return fmt.Errorf("command failed: %w", err)
// 	}

// 	return nil
// }

// func writeLineToFile(filePath, line string) error {
// 	// Open the file in append mode, create it if it doesn't exist
// 	file, err := os.OpenFile(filePath, os.O_APPEND|os.O_CREATE|os.O_WRONLY, 0644)
// 	if err != nil {
// 		return fmt.Errorf("failed to open file: %w", err)
// 	}
// 	defer file.Close()

// 	// Write the line to the file with a newline character
// 	if _, err := file.WriteString(line + "\n"); err != nil {
// 		return fmt.Errorf("failed to write line to file: %w", err)
// 	}

// 	return nil
// }


// resolveBin returns override if non-empty, otherwise the bare command
// name, letting exec.LookPath resolve it from PATH at run time.
func resolveBin(override, fallback string) string {
	if override != "" {
		return override
	}
	return fallback
}

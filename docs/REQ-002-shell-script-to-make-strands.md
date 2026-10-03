# Synthetic Strand Generator Specification (REQ-002)

[← Back to Documentation Index](README.md)

## Objective

Create a shell script (`scripts/make-strands.sh`) that, given a text file,
extracts independent random substrings ("strands") at random lengths and
locations. This simulates the process of shotgun DNA sequencing, creating
synthetic "reads" from a complete reference sequence.

## Inputs

The script accepts the following long-form command-line flags for robust
parsing:

- `-i, --input FILE`: Path to the input text file.
- `-m, --min INT`: Minimum length of a strand.
- `-M, --max INT`: Maximum length of a strand.
- `-n, --count INT`: Number of strands to produce.

The text file must contain a contiguous block of text. The script strips all
newlines before processing to ensure the text is perfectly contiguous.

## Outputs

- Output is written to standard output (`stdout`), allowing simple redirection
  to a file.
- The output contains the extracted strands, formatted with one line for each
  strand.

## Architectural Decisions

### ADR-1: Independent Overlapping Substrings

- **Decision:** Strands are sampled independently and may overlap.
- **Rationale:** This aligns perfectly with the "de novo sequence assembly"
  context of the project, simulating random sequencing fragments rather than
  strictly partitioning the text.

### ADR-2: AWK for Text Extraction

- **Decision:** Use an `awk` block within the Bash script for random sampling.
- **Rationale:** Avoids Bash's `$RANDOM` limit (which caps at 32,767) and
  provides significantly faster text manipulation for large input sequences.

### ADR-3: Newline Stripping

- **Decision:** Strip all newlines from the input text before processing.
- **Rationale:** Guarantees that extracted strands do not contain newlines,
  which would otherwise break the "one line per strand" output format
  requirement.

## See Also

- [Documentation Index](README.md)
- [Assembly Dynamics and Parameter Heuristics](heuristics.md)
- [Assembler Specification (REQ-001)](REQ-001-functional-greedy-overlap-assembler.md)

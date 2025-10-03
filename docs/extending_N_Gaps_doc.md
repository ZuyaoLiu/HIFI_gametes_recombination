# FASTA Gap Expander

Normalize gaps in a genome assembly FASTA by expanding any run of `N` to a fixed length. This helps keep gap representation consistent for downstream analyses and data submission.

## Features
- Works with single- or multi-sequence FASTA files.
- Converts all sequences to uppercase.
- Replaces any contiguous run of `N` (`N+`) with exactly `gapsize` Ns (default: 30,000).
- Wraps output lines to a fixed width (default: 60 characters).
- Optional explicit output path.

## Usage

    python expand_gaps.py input.fasta [--gapsize 30000] [--line-width 60] [-o output.fasta]

### Required arguments
- **input_fa**  
  Path to the input assembly FASTA file.

### Optional arguments
- **--gapsize** (int, default: 30000)  
  Target size to which any run of `N`s is expanded.

- **--line-width** (int, default: 60)  
  Line width for wrapping FASTA output. (Hidden in help, supported by the script.)

- **-o, --output** (Path, optional)  
  Output FASTA path. If omitted, the script may write to stdout or a default location depending on your implementation.

## Examples

    # Expand all gaps in assembly.fasta to 30k Ns and write to expanded.fasta
    python expand_gaps.py assembly.fasta -o expanded.fasta

    # Expand gaps to 10k Ns instead
    python expand_gaps.py assembly.fasta --gapsize 10000 -o expanded_10k.fasta



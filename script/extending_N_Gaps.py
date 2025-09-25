#!/usr/bin/env python3
import argparse
import re
from pathlib import Path
from typing import Iterator


class FastaGapExpander:
    """
    Expand runs of N in a FASTA assembly to a fixed gap size.
    - Reads a single FASTA file (multi-sequence supported).
    - Converts sequences to uppercase.
    - Replaces any contiguous 'N+' with 'N' * gapsize.
    - Wraps output lines to fixed width, default 60.
    """

    def __init__(self, input_fa: Path, gapsize: int = 30000, line_width: int = 60) -> None:
        self.input_fa = Path(input_fa)
        self.gapsize = int(gapsize)
        self.line_width = int(line_width)
        self._n_pattern = re.compile(r"N+")

    def _format_output_path(self) -> Path:
        """
        Build an output filename next to the input.
        Examples:
          gapsize=30000 -> <input>.N30k.fasta
          gapsize=12345 -> <input>.N12345.fasta
        """
        stem = self.input_fa.name
        suffix = ""
        if self.gapsize % 1000 == 0:
            suffix = f".N{self.gapsize // 1000}k.fasta"
        else:
            suffix = f".N{self.gapsize}.fasta"
        return self.input_fa.with_name(stem + suffix)

    def _records(self) -> Iterator[tuple[str, str]]:
        """
        Stream FASTA records as (header, sequence) tuples.
        Sequence is returned as a single uppercase string without whitespace.
        """
        name = None
        seq_chunks = []
        with self.input_fa.open() as fh:
            for line in fh:
                if line.startswith(">"):
                    if name is not None:
                        yield name, "".join(seq_chunks).upper()
                    name = line[1:].strip()
                    seq_chunks = []
                else:
                    seq_chunks.append(line.strip())
            if name is not None:
                yield name, "".join(seq_chunks).upper()

    def _expand_ns(self, s: str) -> str:
        """Replace any run of N with a fixed number of Ns."""
        if not s:
            return s
        return self._n_pattern.sub("N" * self.gapsize, s)

    def _write_wrapped(self, fh, text: str) -> None:
        """Write text wrapped to the configured line width."""
        for i in range(0, len(text), self.line_width):
            fh.write(text[i : i + self.line_width] + "\n")

    def run(self, output_fa: Path | None = None) -> Path:
        """
        Execute the expansion and write to output_fa (or auto-named).
        Returns the output path.
        """
        outpath = Path(output_fa) if output_fa else self._format_output_path()
        with outpath.open("w") as out:
            for name, seq in self._records():
                expanded = self._expand_ns(seq)
                out.write(f">{name}\n")
                self._write_wrapped(out, expanded)
        return outpath


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Expand N gaps in a FASTA assembly to a fixed size (default: 30k)."
    )
    parser.add_argument("input_fa", type=Path, help="Input assembly FASTA file")
    parser.add_argument("--gapsize", type=int, default=30000, help="Target gap size for any N-run (default: 30000)")
    # Hidden option for future flexibility; kept simple per requirements.
    parser.add_argument("--line-width",type=int,default=60,help=argparse.SUPPRESS,)
    # Optional explicit output path (not required by spec, but handy).
    parser.add_argument("-o","--output",type=Path,default=None,help=argparse.SUPPRESS,)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    expander = FastaGapExpander(args.input_fa, args.gapsize, args.line_width)
    outpath = expander.run(args.output)
    print(f"Wrote: {outpath}")


if __name__ == "__main__":
    main()


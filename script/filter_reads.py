#!/usr/bin/env python3
"""
Filter read names that map consistently in two BAMs.

A read is kept if:
  - it appears in BOTH BAMs, and
  - the reference (contig) names are identical, and
  - the strands (orientation) are identical.

By default, only primary alignments are considered (secondary/supplementary are ignored).
"""

from __future__ import annotations
import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Tuple, Iterable
import sys

import pysam


@dataclass(frozen=True)
class MapInfo:
    """Minimal alignment fingerprint used for consistency checking."""
    rname: str
    reverse: bool


class ConsistentReadFilter:
    """Load two BAMs, collect per-read mapping info, and keep reads consistent in both."""

    def __init__(self, bam1: Path, bam2: Path, primary_only: bool = True, min_mapq: int = 0):
        self.bam1 = Path(bam1)
        self.bam2 = Path(bam2)
        self.primary_only = primary_only
        self.min_mapq = int(min_mapq)

    @staticmethod
    def _is_primary(aln: pysam.AlignedSegment) -> bool:
        """Return True if alignment is primary (not secondary or supplementary)."""
        return not (aln.is_secondary or aln.is_supplementary)

    def _collect(self, bam_path: Path) -> Dict[str, MapInfo]:
        """
        Build a dict: read_name -> MapInfo for a BAM.
        If multiple records share the same read name, the last primary record wins.
        """
        info: Dict[str, MapInfo] = {}
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            for aln in bam:
                if aln.is_unmapped:
                    continue
                if self.primary_only and not self._is_primary(aln):
                    continue
                if aln.mapping_quality < self.min_mapq:
                    continue
                rname = bam.get_reference_name(aln.reference_id)
                info[aln.query_name] = MapInfo(rname=rname, reverse=aln.is_reverse)
        return info

    def run(self) -> Tuple[Iterable[str], int]:
        """
        Execute filtering.
        Returns (iterator_of_names, count).
        """
        d1 = self._collect(self.bam1)
        d2 = self._collect(self.bam2)

        keep = []
        # Intersection of read names present in both maps
        for q in (d1.keys() & d2.keys()):
            if d1[q] == d2[q]:
                keep.append(q)

        return keep, len(keep)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Keep read names that map to the same contig and strand in two BAMs."
    )
    p.add_argument("bam1", type=Path, help="First BAM (e.g., hap1.bam)")
    p.add_argument("bam2", type=Path, help="Second BAM (e.g., hap2.bam)")
    p.add_argument(
        "-o", "--output",
        type=Path,
        default=Path("kept.readnames.txt"),
        help="Output file for kept read names (default: kept.readnames.txt)"
    )
    p.add_argument(
        "--no-primary-only",
        dest="primary_only",
        action="store_false",
        help="Consider secondary/supplementary alignments as well (default: primary only)."
    )
    p.add_argument(
        "--min-mapq",
        type=int,
        default=0,
        help="Minimum mapping quality to consider (default: 0)."
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()

    # Basic input validation
    for bam in (args.bam1, args.bam2):
        if not bam.exists():
            print(f"ERROR: BAM not found: {bam}", file=sys.stderr)
            sys.exit(2)
        if bam.suffix not in (".bam", ".cram"):
            print(f"WARNING: {bam} does not look like a BAM/CRAM file.", file=sys.stderr)

    try:
        filt = ConsistentReadFilter(
            args.bam1, args.bam2,
            primary_only=args.primary_only,
            min_mapq=args.min_mapq
        )
        names, n = filt.run()
    except Exception as e:
        print(f"ERROR: failed during filtering: {e}", file=sys.stderr)
        sys.exit(3)

    try:
        with open(args.output, "w") as w:
            for q in names:
                w.write(q + "\n")
    except Exception as e:
        print(f"ERROR: cannot write output file {args.output}: {e}", file=sys.stderr)
        sys.exit(4)

    print(f"Kept {n} reads -> {args.output}")


if __name__ == "__main__":
    main()


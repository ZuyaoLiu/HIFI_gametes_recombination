#!/usr/bin/env python3
"""
Joint-partition SNP caller with flank checks on both haplotypes.

A read is considered only if it exists in BOTH BAMs. For each shared read:
- Discard if total soft-clip bases > --softclip-max in either hap.
- Iterate per-read positions; consider only '=' (match) or 'X' (mismatch) on BOTH BAMs at the same read index.
- Count flank matches '=' to the left and right; require >= --flank-len on each side per hap.
- Emit TSV rows for positions where one hap is '=' and the other is 'X'.
- Discard the entire read if positions with (X/X) exceed --bothmm-max.

TSV (stdout) columns:
  read_id            Read name (query name) from the FASTQ/BAM.
  read_length        Length of the read sequence in bases.
  SNP_pos_on_read    0-based index of the candidate SNP position on the read axis.
  baseQ              Phred-scaled base quality at this read position.

  hap1_chr           Contig (reference name) where hap1 alignment is placed.
  hap1_refpos        0-based reference coordinate of hap1 alignment at this read position,
                     or -1 if not aligned to a specific reference base (e.g. insertion/softclip).
  hap1_op            CIGAR-derived operation symbol at this read position for hap1:
                     '=' (match), 'X' (mismatch), 'I' (insertion), 'S' (soft-clip).

  hap2_chr           Contig (reference name) where hap2 alignment is placed.
  hap2_refpos        0-based reference coordinate of hap2 alignment at this read position,
                     or -1 if not aligned.
  hap2_op            CIGAR-derived operation symbol for hap2 at this read position.

  which_hap_matches  Label indicating which haplotype agrees with the reference:
                     "hap1" if hap1 is '=', hap2 is 'X'; "hap2" if hap2 is '=', hap1 is 'X'.

  hap1_flankL        Count of consecutive '=' bases immediately to the LEFT of this position
                     on the hap1 read axis.
  hap1_flankR        Count of consecutive '=' bases immediately to the RIGHT of this position
                     on the hap1 read axis.
  hap2_flankL        Same as hap1_flankL but for hap2.
  hap2_flankR        Same as hap1_flankR but for hap2.

  hap1_flank_ok      Boolean flag (True/False). True if both hap1_flankL and hap1_flankR
                     are >= required flank length (default 10).
  hap2_flank_ok      Boolean flag for hap2 flank check.
  both_flank_ok      Boolean flag. True if BOTH haplotypes individually satisfy the flank check.

Metrics (stderr):
  total_shared_reads, softclip_filtered, bothmm_filtered,
  kept_reads, reads_with_candidates, total_candidate_rows
"""

from __future__ import annotations
import argparse
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple, Dict

import pysam


# ------------------------------- utilities ------------------------------------

def total_softclip(aln: pysam.AlignedSegment) -> int:
    """Sum over all 'S' (soft-clip) CIGAR operations."""
    total = 0
    if aln.cigartuples:
        for op, length in aln.cigartuples:
            if op == 4:  # 'S'
                total += length
    return total


def build_read_maps(
    aln: pysam.AlignedSegment,
) -> Tuple[List[str], List[Optional[int]], str]:
    """
    Expand an alignment onto the read axis.

    Returns:
      ops:    list of symbols per read index: '=', 'X', 'I', 'S', '?'
      refpos: list of 0-based reference positions for '='/'X', else None
      refname: reference/contig name
    """
    qlen = aln.query_length or 0
    ops: List[str] = ['?'] * qlen
    refpos: List[Optional[int]] = [None] * qlen
    refname = aln.reference_name or "*"

    q = 0
    r = aln.reference_start if aln.reference_start is not None else 0

    # pysam CIGAR codes: 0:M,1:I,2:D,3:N,4:S,5:H,6:P,7:=,8:X
    for op, length in (aln.cigartuples or []):
        if op == 7:  # '='
            for _ in range(length):
                if q < qlen:
                    ops[q] = '='
                    refpos[q] = r
                q += 1
                r += 1
        elif op == 8:  # 'X'
            for _ in range(length):
                if q < qlen:
                    ops[q] = 'X'
                    refpos[q] = r
                q += 1
                r += 1
        elif op == 0:  # 'M' (fallback if aligner didn't emit =/X)
            for _ in range(length):
                if q < qlen:
                    ops[q] = '='
                    refpos[q] = r
                q += 1
                r += 1
        elif op == 1:  # 'I'
            for _ in range(length):
                if q < qlen:
                    ops[q] = 'I'
                    refpos[q] = None
                q += 1
        elif op == 2:  # 'D' (consumes ref but not read)
            r += length
        elif op == 4:  # 'S'
            for _ in range(length):
                if q < qlen:
                    ops[q] = 'S'
                    refpos[q] = None
                q += 1
        elif op == 3:  # 'N' (ref skip)
            r += length
        # 5:H and 6:P do not consume read or ref in this context

    return ops, refpos, refname


def flank_match_lengths(ops: List[str], i: int) -> Tuple[int, int]:
    """
    Count contiguous '=' to the LEFT and RIGHT of position i on the read axis.
    The base at index i itself is not included.
    """
    # left
    L = 0
    j = i - 1
    while j >= 0 and ops[j] == '=':
        L += 1
        j -= 1
    # right
    R = 0
    n = len(ops)
    j = i + 1
    while j < n and ops[j] == '=':
        R += 1
        j += 1
    return L, R


# ------------------------------ core classes ----------------------------------

@dataclass
class Thresholds:
    softclip_max: int = 200    # discard read if total S > this (either hap)
    bothmm_max: int = 100      # discard read if (# of X/X positions) > this
    flank_len: int = 10        # require >= this many '=' on both sides per hap


class JointPartitionSNPCaller:
    """Compare two BAMs per-read and emit candidate SNP rows where only one hap mismatches."""

    def __init__(self, bam1: Path, bam2: Path, thr: Thresholds):
        self.bam1 = Path(bam1)
        self.bam2 = Path(bam2)
        self.thr = thr

        # metrics
        self.total_shared_reads = 0
        self.softclip_filtered = 0
        self.bothmm_filtered = 0
        self.kept_reads = 0
        self.reads_with_candidates = 0
        self.total_candidate_rows = 0

    def _load_primary_records(self, bam_path: Path) -> Dict[str, pysam.AlignedSegment]:
        """
        Load one record per read name. If multiple records appear for the same read,
        the last one wins (consistent with the original script's behavior).
        """
        records: Dict[str, pysam.AlignedSegment] = {}
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            for rec in bam:
                if rec.is_unmapped:
                    continue
                records[rec.query_name] = rec
        return records

    def _emit_header(self) -> None:
        cols = [
            "read_id", "read_length", "SNP_pos_on_read", "baseQ",
            "hap1_chr", "hap1_refpos", "hap1_op",
            "hap2_chr", "hap2_refpos", "hap2_op",
            "which_hap_matches",
            "hap1_flankL", "hap1_flankR", "hap2_flankL", "hap2_flankR",
            "hap1_flank_ok", "hap2_flank_ok", "both_flank_ok",
        ]
        print("\t".join(cols))

    def _print_metrics(self) -> None:
        lines = [
            f"[metrics] total_shared_reads={self.total_shared_reads}",
            f"[metrics] softclip_filtered={self.softclip_filtered}",
            f"[metrics] bothmm_filtered={self.bothmm_filtered}",
            f"[metrics] kept_reads={self.kept_reads}",
            f"[metrics] reads_with_candidates={self.reads_with_candidates}",
            f"[metrics] total_candidate_rows={self.total_candidate_rows}",
        ]
        sys.stderr.write("\n".join(lines) + "\n")

    def run(self) -> int:
        """Main entry: load BAMs, iterate shared reads, emit TSV rows, print metrics."""
        # input validation
        for p in (self.bam1, self.bam2):
            if not p.exists():
                sys.stderr.write(f"ERROR: BAM not found: {p}\n")
                return 2

        m1 = self._load_primary_records(self.bam1)
        m2 = self._load_primary_records(self.bam2)

        shared = list(m1.keys() & m2.keys())
        self.total_shared_reads = len(shared)

        self._emit_header()

        for qname in shared:
            a1 = m1[qname]
            a2 = m2[qname]

            # soft-clip filter
            if total_softclip(a1) > self.thr.softclip_max or total_softclip(a2) > self.thr.softclip_max:
                self.softclip_filtered += 1
                continue

            ops1, ref1, chr1 = build_read_maps(a1)
            ops2, ref2, chr2 = build_read_maps(a2)

            qlen = min(len(ops1), len(ops2))
            quals = a1.query_qualities
            has_qual = quals is not None

            rows: List[str] = []
            both_mismatch_errors = 0

            for i in range(qlen):
                op1 = ops1[i]
                op2 = ops2[i]

                # skip soft-clipped bases
                if op1 == 'S' or op2 == 'S':
                    continue

                # only consider positions that consume the read base in BOTH alignments
                if (op1 in ('=', 'X')) and (op2 in ('=', 'X')):
                    if op1 == 'X' and op2 == 'X':
                        both_mismatch_errors += 1
                        continue

                    if op1 == '=' and op2 == 'X':
                        which = "hap1"
                    elif op1 == 'X' and op2 == '=':
                        which = "hap2"
                    else:
                        continue  # '=' vs '='

                    # per-hap flank checks
                    h1L, h1R = flank_match_lengths(ops1, i)
                    h2L, h2R = flank_match_lengths(ops2, i)
                    h1_ok = (h1L >= self.thr.flank_len and h1R >= self.thr.flank_len)
                    h2_ok = (h2L >= self.thr.flank_len and h2R >= self.thr.flank_len)
                    both_ok = (h1_ok and h2_ok)

                    baseQ = int(quals[i]) if has_qual else 0
                    hap1_refpos = ref1[i] if ref1[i] is not None else -1
                    hap2_refpos = ref2[i] if ref2[i] is not None else -1

                    rows.append("\t".join(map(str, [
                        qname,
                        qlen,
                        i,
                        baseQ,
                        chr1,
                        hap1_refpos,
                        op1,
                        chr2,
                        hap2_refpos,
                        op2,
                        which,
                        h1L,
                        h1R,
                        h2L,
                        h2R,
                        h1_ok,
                        h2_ok,
                        both_ok,
                    ])))

            # whole-read filter on X/X overload
            if both_mismatch_errors > self.thr.bothmm_max:
                self.bothmm_filtered += 1
                continue

            self.kept_reads += 1

            if rows:
                self.reads_with_candidates += 1
                self.total_candidate_rows += len(rows)
                print("\n".join(rows))

        self._print_metrics()
        return 0


# ---------------------------------- CLI ---------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Joint-partition SNP caller with per-read flank checks on two BAMs."
    )
    p.add_argument("bam1", type=Path, help="First BAM (e.g. hap1.bam)")
    p.add_argument("bam2", type=Path, help="Second BAM (e.g. hap2.bam)")
    p.add_argument("--softclip-max", type=int, default=200,
                   help="Max total soft-clip bases allowed per read in each BAM (default: 200)")
    p.add_argument("--bothmm-max", type=int, default=100,
                   help="Max allowed X/X positions per read before discarding the read (default: 100)")
    p.add_argument("--flank-len", type=int, default=10,
                   help="Required '=' flank length on BOTH sides per hap (default: 10)")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    thr = Thresholds(
        softclip_max=args.softclip_max,
        bothmm_max=args.bothmm_max,
        flank_len=args.flank_len,
    )
    caller = JointPartitionSNPCaller(args.bam1, args.bam2, thr)
    rc = caller.run()
    sys.exit(rc)


if __name__ == "__main__":
    main()


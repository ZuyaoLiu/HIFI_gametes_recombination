#!/usr/bin/env python3
"""
Map SNP positions on reads (given by positions along the read sequence) to reference coordinates,
using an aligned BAM (e.g., PacBio HiFi reads aligned with minimap2).

Input TSV schema (tab-delimited; header required):
    read_id SNP1_pos_on_read SNP2_pos_on_read CO1_hap1_chr CO1_hap1_refpos CO2_hap1_chr CO2_hap1_refpos CO1_hap2_chr CO1_hap2_refpos CO2_hap2_chr CO2_hap2_refpos Event

Notes:
- SNP*_pos_on_read are interpreted as 1-based coordinates along the read sequence (as shown by many variant callers and common tables).
- NA is allowed in any column; NA in SNP positions will be skipped.
- The script will map SNP1 and SNP2 (when present) to reference (chrom, pos, strand) using the read's *primary* alignment by default.
- If the read has multiple alignments, the script uses: primary & not supplementary preferred; otherwise highest MAPQ; ties broken by longest aligned ref span.
- Coordinates on the reference are output as 1-based positions (SAM convention for POS).

Output TSV adds columns (after the original ones):
    aln_rname aln_pos aln_end mapq cigar strand is_secondary is_supplementary
    SNP1_ref_chrom SNP1_ref_pos SNP1_status
    SNP2_ref_chrom SNP2_ref_pos SNP2_status

Where status ∈ {mapped, not_aligned, clipped, insertion_gap, deletion_gap, out_of_range, NA}.
"""

import sys
import argparse
import pandas as pd
import pysam

def pick_best_alignment(alignments):
    """Prefer primary (not secondary, not supplementary). Else highest MAPQ; tie-breaker: longest ref span."""
    if not alignments:
        return None
    primaries = [a for a in alignments if not a.is_secondary and not a.is_supplementary]
    cands = primaries if primaries else alignments
    # Highest MAPQ
    max_mapq = max(a.mapping_quality for a in cands)
    cands = [a for a in cands if a.mapping_quality == max_mapq]
    if len(cands) == 1:
        return cands[0]
    # Tie-breaker: longest aligned reference span
    def ref_span(a):
        return a.reference_end - a.reference_start if a.reference_start is not None and a.reference_end is not None else -1
    cands.sort(key=ref_span, reverse=True)
    return cands[0]
def build_read_to_ref_map(aln):
    """
    Build a mapping from read positions (0-based) to reference positions (1-based).
    Uses get_aligned_pairs(matches_only=False), which returns tuples (read_pos, ref_pos).
      - read_pos: 0-based coordinate on the read, or None if it's a deletion in the read
      - ref_pos: 0-based coordinate on the reference, or None if it's an insertion/clip
    """
    mapping = {}
    for read_pos, ref_pos in aln.get_aligned_pairs(matches_only=False):
        if read_pos is None:
            # Deletion in the read relative to reference (no base in read to map)
            continue
        if ref_pos is None:
            # Insertion or soft-clip relative to reference: no ref coordinate
            mapping[read_pos] = None
        else:
            # Convert ref_pos from 0-based (pysam) to 1-based (standard SAM convention)
            mapping[read_pos] = ref_pos + 1
    return mapping
def map_read_position(aln, read_pos_1based):
    """
    Map a 1-based position along the read to a reference coordinate.
    Returns: (chrom, ref_pos, status, strand)

    Status values:
      - "mapped": the read base maps directly to a reference base
      - "not_aligned": the read is not aligned at all, or the position is outside alignment
      - "out_of_range": the SNP position is outside the read length
      - "nearest_mapped_from_insertion_or_clipped(offset=±N)":
            the read base falls inside an insertion or soft-clip,
            so the closest mapped neighbor base was used instead,
            with offset = distance in read bases (negative = left, positive = right)
      - "NA": input position was NA
    """
    if read_pos_1based is None:
        return (None, None, "NA", "+")
    if aln.reference_id is None:
        return (None, None, "not_aligned", "+")

    strand = "-" if aln.is_reverse else "+"
    read_index = read_pos_1based - 1

    # Get read length
    qlen = aln.query_length
    if qlen is None:
        qseq = aln.query_sequence or ""
        qlen = len(qseq)

    if read_index < 0 or (qlen is not None and read_index >= qlen):
        return (None, None, "out_of_range", strand)

    r2r = build_read_to_ref_map(aln)

    # Direct mapping
    ref_pos = r2r.get(read_index, None)
    if ref_pos is not None:
        return (aln.reference_name, ref_pos, "mapped", strand)

    # If read_index is absent from r2r: not aligned at all
    if read_index not in r2r:
        return (None, None, "not_aligned", strand)

    # Insertion or soft-clip: search nearest mapped neighbor
    nearest_pos = None
    nearest_offset = None
    for k in range(1, qlen):  # expand search outward
        left = read_index - k
        right = read_index + k
        candidates = []

        if left >= 0 and left in r2r and r2r[left] is not None:
            candidates.append(("L", left, r2r[left], -k))
        if right < qlen and right in r2r and r2r[right] is not None:
            candidates.append(("R", right, r2r[right], +k))

        if candidates:
            # Take the first found: left first, then right (tie broken by order)
            side, idx, refp, off = candidates[0]
            nearest_pos = refp
            nearest_offset = off
            break

    if nearest_pos is not None:
        status = f"nearest_mapped_from_insertion_or_clipped(offset={nearest_offset:+d})"
        return (aln.reference_name, nearest_pos, status, strand)

    # No mapped neighbors found (very rare case)
    return (None, None, "not_aligned", strand)

def parse_pos(x):
    if isinstance(x, (int, float)):
        if pd.isna(x):
            return None
        # allow floats like 2831566.5 -> round? Typically SNP positions should be integers. We'll round to nearest.
        return int(round(x))
    if isinstance(x, str):
        xs = x.strip()
        if xs == "" or xs.upper() == "NA":
            return None
        try:
            if "." in xs:
                return int(round(float(xs)))
            return int(xs)
        except ValueError:
            return None
    return None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True, help="Aligned HiFi reads BAM (indexed).")
    ap.add_argument("--tsv", required=True, help="Input TSV with read_id and SNP*_pos_on_read.")
    ap.add_argument("--out", required=True, help="Output TSV path.")
    args = ap.parse_args()

    bam = pysam.AlignmentFile(args.bam, "rb")
    df = pd.read_csv(args.tsv, sep="\t", dtype=str, keep_default_na=False)
    # Normalize columns
    expected_cols = ["read_id","SNP1_pos_on_read","SNP2_pos_on_read",
                     "CO1_hap1_chr","CO1_hap1_refpos","CO2_hap1_chr","CO2_hap1_refpos",
                     "CO1_hap2_chr","CO1_hap2_refpos","CO2_hap2_chr","CO2_hap2_refpos","Event"]
    missing = [c for c in expected_cols if c not in df.columns]
    if missing:
        print(f"WARNING: missing columns in TSV: {missing}", file=sys.stderr)

    # Prepare outputs
    out_rows = []
    # Build index: read_id -> list of alignments
    # We'll fetch by read name using the BAM index; not all BAMs support name index, so we scan once.
    # For reasonably sized BAMs this is fine; for huge BAMs, recommend name-sorted BAM with an auxiliary index.
    aln_by_name = {}
    for aln in bam.fetch(until_eof=True):
        rn = aln.query_name
        if rn is None:
            continue
        aln_by_name.setdefault(rn, []).append(aln)
    for _, row in df.iterrows():
        read_id = row.get("read_id", "").strip()
        snp1_read = parse_pos(row.get("SNP1_pos_on_read", ""))
        snp2_read = parse_pos(row.get("SNP2_pos_on_read", ""))
        aligns = aln_by_name.get(read_id, [])
        best = pick_best_alignment(aligns)
        if best is None:
            aln_rname = aln_pos = aln_end = mapq = cigar = strand = None
            is_sec = is_sup = None
            s1c, s1p, s1s, st = (None, None, "not_aligned", "+")
            s2c, s2p, s2s, _ = (None, None, "not_aligned", "+")
        else:
            aln_rname = best.reference_name
            aln_pos = (best.reference_start + 1) if best.reference_start is not None else None
            aln_end = best.reference_end if best.reference_end is not None else None
            mapq = best.mapping_quality
            cigar = best.cigarstring
            strand = "-" if best.is_reverse else "+"
            is_sec = best.is_secondary
            is_sup = best.is_supplementary
            s1c, s1p, s1s, _ = map_read_position(best, snp1_read)
            s2c, s2p, s2s, _ = map_read_position(best, snp2_read)

        new = {
            "read_id": row.get("read_id", ""),
            "SNP1_pos_on_read": row.get("SNP1_pos_on_read", ""),
            "SNP2_pos_on_read": row.get("SNP2_pos_on_read", ""),
            "Event": row.get("Event", ""),
        
            "aln_rname": aln_rname,
            "aln_pos": aln_pos,
            "aln_end": aln_end,
            "mapq": mapq,
        
            "SNP1_ref_chrom": s1c,
            "SNP1_ref_pos": s1p,
            "SNP1_status": s1s,
        
            "SNP2_ref_chrom": s2c,
            "SNP2_ref_pos": s2p,
            "SNP2_status": s2s,
        }
        out_rows.append(new)
    out_df = pd.DataFrame(out_rows)
    out_df.to_csv(args.out, sep="\t", index=False)

if __name__ == "__main__":
    main()

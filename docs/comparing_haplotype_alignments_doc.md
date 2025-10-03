# Joint-Partition SNP Caller with Flank Checks

This tool performs **joint SNP calling** from two haplotype-aligned BAM files. It compares per-read alignments of the same sequencing read across both haplotypes and identifies candidate SNPs where one haplotype matches the reference and the other shows a mismatch.  

The method enforces quality checks such as soft-clipping limits, flanking match requirements, and mismatch thresholds to minimize false positives.

---

## Overview

- Only reads present in **both BAMs** are considered.  
- Reads are filtered if:
  - The total number of soft-clipped bases exceeds `--softclip-max` in **either** BAM.  
  - The number of `(X/X)` mismatches (both haplotypes disagree with the reference at the same position) exceeds `--bothmm-max`.  
- Within each read:
  - Positions are considered only if both haplotypes have a valid alignment operation (`=` or `X`) at the same read index.  
  - A candidate SNP is reported if one haplotype has `=` (match) while the other has `X` (mismatch).  
  - For each candidate SNP, the flanking bases are checked: each haplotype must have at least `--flank-len` consecutive `=` matches both to the left and right of the SNP position.  

---

## Output

### Main Output (TSV to **stdout**)

Each row corresponds to one candidate SNP position. Columns:

- **read_id** – Read name (query name).  
- **read_length** – Length of the read sequence in bases.  
- **SNP_pos_on_read** – 0-based index of the SNP position on the read axis.  
- **baseQ** – Phred-scaled base quality at the SNP position.  

**Haplotype 1 information**:
- **hap1_chr** – Contig name where the read is aligned.  
- **hap1_refpos** – 0-based reference coordinate at this read position. 
- **hap1_op** – Operation symbol from the CIGAR string at this position:  
  - `=` (match), `X` (mismatch)

**Haplotype 2 information**:
- **hap2_chr** – Contig name.  
- **hap2_refpos** – 0-based reference coordinate.  
- **hap2_op** – Operation symbol from the CIGAR string at this position:  
  - `=` (match), `X` (mismatch)  

**Comparison and flank metrics**:
- **which_hap_matches** – Indicates which haplotype agrees with the reference:  
  - `"hap1"` if hap1 is `=`, hap2 is `X`.  
  - `"hap2"` if hap2 is `=`, hap1 is `X`.  
- **hap1_flankL / hap1_flankR** – Number of consecutive `=` matches immediately to the left/right of this position on hap1.  
- **hap2_flankL / hap2_flankR** – Same metrics for hap2.  
- **hap1_flank_ok** – Boolean (`True`/`False`), `True` if both hap1_flankL and hap1_flankR ≥ required flank length.  
- **hap2_flank_ok** – Boolean flag for hap2 flank check.  
- **both_flank_ok** – `True` only if both haplotypes individually satisfy the flank check.  

---

### Metrics (written to **stderr**)

Summary counts printed to stderr:

- **total_shared_reads** – Number of reads found in both BAMs.  
- **softclip_filtered** – Reads discarded due to exceeding the soft-clip threshold.  
- **bothmm_filtered** – Reads discarded due to too many `(X/X)` mismatches.  
- **kept_reads** – Reads that passed all filters.  
- **reads_with_candidates** – Reads that produced at least one candidate SNP row.  
- **total_candidate_rows** – Total number of SNP candidate positions output.  

---

## Usage

    python joint_snp_caller.py hap1.bam hap2.bam [options]

### Required arguments
- **bam1**  
  First BAM file (e.g., `hap1.bam`).  

- **bam2**  
  Second BAM file (e.g., `hap2.bam`).  

### Optional arguments
- **--softclip-max** (int, default: 200)  
  Maximum total number of soft-clipped bases allowed per read in **each** BAM.  

- **--bothmm-max** (int, default: 100)  
  Maximum number of `(X/X)` mismatch positions allowed per read. Reads exceeding this are discarded.  

- **--flank-len** (int, default: 10)  
  Required number of consecutive `=` matches flanking both sides of a candidate SNP, for **each haplotype**.  

---

## Examples

1. **Run with default parameters (softclip max 200, flank length 10):**

    ```bash
    python joint_snp_caller.py hap1.bam hap2.bam > candidates.tsv 2> metrics.log
    ```

2. **Increase flank stringency (20 matches required on each side):**

    ```bash
    python joint_snp_caller.py hap1.bam hap2.bam --flank-len 20 > candidates.tsv
    ```

3. **Filter out reads with more than 50 X/X mismatches:**

    ```bash
    python joint_snp_caller.py hap1.bam hap2.bam --bothmm-max 50
    ```

---



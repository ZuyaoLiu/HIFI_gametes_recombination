#  Mapping Read-Based SNP Positions to Reference Coordinates

## Overview
This script maps SNP positions recorded **relative to the read sequence** (positions along the read, e.g. 1–N) to their corresponding **reference genome coordinates**, using an aligned BAM file (e.g. PacBio HiFi reads aligned with minimap2).

It is designed to post-process candidate recombination events or SNP transitions previously called on read sequences, and add genomic coordinates plus alignment context.

---

## Input

### 1. BAM file
- Contains alignments of reads to the reference genome.
- Must be coordinate- or name-sorted and **indexed**.
- Script scans through the BAM to build an index of `read_id → alignments`.

### 2. TSV file
- Tab-delimited with at least these columns (header required):
  - `read_id`  
  - `SNP1_pos_on_read`, `SNP2_pos_on_read`: 1-based SNP positions along the read sequence  
  - `CO1_hap1_chr`, `CO1_hap1_refpos`, … `CO2_hap2_refpos` (other columns from upstream analysis; optional)  
  - `Event`: classification (e.g., CO/NCO/Complex)  

Values may be `NA`. SNP positions may be integer or float (floats are rounded).

---

## Output
A new TSV file with the original input columns **plus additional alignment- and mapping-related columns**:

- **Alignment summary**  
  - `aln_rname`: reference contig name  
  - `aln_pos`: alignment start position (1-based)  
  - `aln_end`: alignment end position  
  - `mapq`: mapping quality  
  - `strand`: `+` or `-`  
  - `is_secondary`, `is_supplementary`: flags for alignment type  

- **SNP 1 mapping**  
  - `SNP1_ref_chrom`, `SNP1_ref_pos`: mapped genomic coordinate  
  - `SNP1_status`: mapping status  

- **SNP 2 mapping**  
  - `SNP2_ref_chrom`, `SNP2_ref_pos`  
  - `SNP2_status`  

### Status values
- `mapped` : read base maps directly to reference  
- `not_aligned` : read not aligned or SNP position outside alignment  
- `out_of_range` : SNP position exceeds read length  
- `nearest_mapped_from_insertion_or_clipped(offset=±N)` : SNP fell inside insertion/clip; mapped to nearest aligned neighbor, with offset in read bases  
- `NA` : input SNP position was missing  

---

## Workflow Details

### 1. Pick best alignment per read
- If multiple alignments exist, choose:
  1. Primary alignment (not secondary, not supplementary), else all alignments  
  2. Among candidates, highest MAPQ  
  3. If tie: longest aligned reference span  

### 2. Build read→reference map
- Using `aln.get_aligned_pairs(matches_only=False)`  
- Constructs dictionary of read positions (0-based) to reference positions (1-based)  
- Handles deletions, insertions, and clipping  

### 3. Map SNP positions
- For each SNP position (1-based along read):
  - If mapped directly → `mapped`  
  - If insertion/clip → search outward for nearest aligned base → `nearest_mapped_from_insertion_or_clipped`  
  - If no mapping possible → `not_aligned`  

### 4. Handle edge cases
- `out_of_range` if SNP position <1 or > read length  
- Floats or non-integers are rounded to nearest integer  

### 5. Append results
- Original row is preserved with new alignment/mapping columns appended  
- Output is written as TSV  

---
### Complex Event Logic

When a read shows **more than two haplotype transitions**, it is classified as a **Complex** event.  
In this case:  

- The script only records the **first two detected transitions** (similar to the NCO case).  
- For each of these two transitions, the output includes:
  - `SNP1_pos_on_read` and `SNP2_pos_on_read`: read-level positions of the SNPs where the first and second haplotype switches occur.  
  - Midpoints of the flanking haplotype reference coordinates (`hap1_refpos`, `hap2_refpos`) for both transitions.  

Thus, although the read may contain **three or more switches**, only the first two are output, and the event type is labeled as `"Complex"`.  
This preserves partial coordinate information but signals that the underlying haplotype pattern cannot be explained by a simple CO or NCO model.  

---

## Usage
```bash
python map_snp_to_ref.py --bam reads.bam --tsv events.tsv --out events_with_ref.tsv

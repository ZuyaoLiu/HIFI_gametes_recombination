# Haplotype-Based Candidate Recombinant Read Filtering

## Overview
This script processes read-level SNP–haplotype matching data to identify **candidate recombinant reads** while filtering out spurious signals. It is designed to work on a tab-delimited file where each row corresponds to a SNP observed on a sequencing read, with haplotype-matching information.  

The pipeline:  
1. Assigns haplotype ratios to each read.  
2. Identifies “pure” reads that consistently support one haplotype.  
3. Filters SNP sites to keep only those supported by sufficient pure reads.  
4. Re-assigns haplotype ratios on the filtered data.  
5. Flags candidate reads showing haplotype transitions.  
6. Applies two filters to candidate reads:  
   - **Uniqueness**: transition pairs must be unique across reads.  
   - **Coverage**: each interval around a transition must be supported by sufficient pure reads.  
7. Outputs filtered candidate reads, dropped reads (with reasons), and transition pair information.  

---

## Input
- **Tab-delimited file (`input_file.tsv`)** with at least the following columns:
  - `read_id`: identifier of each sequencing read  
  - `SNP_pos_on_read`: SNP position along the read  
  - `which_hap_matches`: `"hap1"` or `"hap2"`, indicating haplotype match  
  - `hap1_chr`, `hap1_refpos`: coordinates of the SNP on haplotype 1  
  - `hap2_chr`, `hap2_refpos`: coordinates of the SNP on haplotype 2  

Additional columns are preserved and used where relevant.  

---

## Output
For a given prefix (e.g. `results`), the script produces:  

- `results_final_filt.tsv`  
  Candidate recombinant reads that passed all filters (per-SNP rows).  

- `results_dropped.tsv`  
  Information about reads dropped due to recurring transitions or insufficient pure read coverage.  

- `results_trans_pairs.tsv`  
  Transition pair table, showing haplotype switch events detected per read, annotated with coverage and recurrence statistics.  

---

## Workflow Details

### 1. Haplotype assignment (`haplotype_assignment`)
- Groups SNPs by read.  
- Counts number of SNPs matching `hap1` vs `hap2`.  
- Calculates per-read haplotype ratios (`hap1_ratio`, `hap2_ratio`).  

### 2. Pure read selection (`filter_reads_95consistency`)
- Defines **pure reads** as those with ≥95% consistency for one haplotype.  
- Returns lists of read IDs for hap1-pure and hap2-pure.  

### 3. SNP filtering (`filter_SNP_using_pure_reads`)
- Retains SNPs only if both haplotypes are supported by ≥3 pure reads.  
- Produces a filtered dataframe restricted to confident SNP positions.  

### 4. Recalculate haplotype ratios  
- Runs `haplotype_assignment` again on the filtered SNP set.  

### 5. Candidate read identification (`identify_candidate_reads`)
- Candidate reads must show mixed haplotype support:  
  - hap1_ratio between 0.05 and 0.95  
  - hap2_ratio between 0.05 and 0.95  
- Marks these reads in the dataframe (`is_candidate = True`).  

### 6. Candidate filtering 
- **Transition pair construction**:  
  - Detects adjacent SNPs within a read where haplotype switches occur.  
  - Records the transition interval coordinates on both haplotypes.  
- **Rule A – Uniqueness**: transition pairs recurring across different reads are flagged.  
- **Rule B – Coverage**: for each transition interval, requires at least `min_pure_cov` pure reads covering both haplotypes (default: 3).  
- Reads with recurring pairs or insufficient coverage are dropped.  

### 7. Output results
- Keeps only candidate reads that pass both rules.  
- Saves final filtered candidates, dropped read info, and transition-pair audit table.  

---

## Usage
```bash
python script.py input_file.tsv output_prefix

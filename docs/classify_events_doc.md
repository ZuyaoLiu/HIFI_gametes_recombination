# Read-Level Crossover/Non-Crossover Classification

## Overview
This script analyzes sequencing reads aligned to phased haplotypes and attempts to classify each read into one of several categories of **recombination events** based on observed haplotype transitions between adjacent SNPs:  

- **CO (Crossover)**: a single clean haplotype switch in the interior of the read.  
- **NCO (Non-Crossover / Gene Conversion)**: two clear haplotype switches (entry and exit), producing a short conversion tract.  
- **Ambiguous**: a single switch but located too close to the read ends (uncertain classification).  
- **Complex**: more than two haplotype switches, suggesting complicated or noisy events.  


---

## Input
- **Tab-delimited file** with at least these columns:
  - `read_id` : identifier of sequencing read  
  - `SNP_pos_on_read` : SNP position along the read (integer)  
  - `which_hap_matches` : `"hap1"` or `"hap2"` indicating haplotype match at that SNP  
  - `hap1_chr`, `hap1_refpos` : haplotype 1 genomic coordinates of the SNP  
  - `hap2_chr`, `hap2_refpos` : haplotype 2 genomic coordinates of the SNP  

The file should contain multiple SNP entries per read.  

---

## Output
- `rec_events.tsv`  
  Table of per-read recombination classifications. Columns:

  - `read_id`: read identifier  
  - `SNP1_pos_on_read`, `SNP2_pos_on_read`: read-level SNP positions flanking the transition(s)  
  - `CO1_hap1_chr`, `CO1_hap1_refpos`: hap1 coordinates near first switch  
  - `CO2_hap1_chr`, `CO2_hap1_refpos`: hap1 coordinates near second switch (if present)  
  - `CO1_hap2_chr`, `CO1_hap2_refpos`: hap2 coordinates near first switch  
  - `CO2_hap2_chr`, `CO2_hap2_refpos`: hap2 coordinates near second switch (if present)  
  - `Event`: event type (`CO`, `NCO`, `Ambiguous`, `Complex`)  

---

## Workflow Details

### 1. Classify SNP haplotype per read
- For each read, the script collects the sequence of haplotype labels (`hap1`=1, `hap2`=2) across SNPs.  

### 2. Detect transitions
- A **transition** is defined as two adjacent SNPs with different haplotype matches (`hap1` → `hap2` or vice versa).  
- All transitions are recorded; their count determines classification.  

### 3. Classification rules
- **0 transitions**: read is pure haplotype (ignored in output).  
- **1 transition**:
  - If it occurs **not at the edges** of the read (SNP index ≥1 and ≤ n-2): → classified as **CO**.  
  - If it is **too close to the edges**: → classified as **Ambiguous**.  
- **2 transitions**: interpreted as **NCO (non-crossover / conversion tract)**.  
- **>2 transitions**: interpreted as **Complex**.  

### 4. Event coordinate reporting
- For each transition, the script reports the midpoints of the flanking SNPs’ haplotype reference positions (`hap1_refpos`, `hap2_refpos`).  
- These serve as approximate genomic breakpoints for the recombination tract.  

### 5. Results aggregation
- One row per read is returned, including SNP positions, haplotype coordinates, and the assigned `Event`.  
- Output is sorted by `CO1_hap1_chr` and `CO1_hap1_refpos`.  

---

## Usage
```bash
python script.py input_file.tsv

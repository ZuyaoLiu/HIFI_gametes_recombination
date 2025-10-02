# Identifying Crossovers and Gene Conversions from Gamete HiFi Sequencing (PacBio Revio)

This workflow is designed to **detect crossovers and gene conversion events** from PacBio HiFi sequencing (Revio platform) by comparing phased haploid assemblies and read alignments. It applies stringent filtering on read mappings and SNP calls to ensure accurate detection of recombination signals.  

While originally developed for **gamete sequencing** (e.g., sperm), the pipeline is also generalizable to **bulk sequencing data** (such as cancer tumor samples), where distinguishing haplotype-specific variants and identifying recombination-like events can reveal insights into **genomic instability and clonal evolution**.  

---

## Key Features

- **Scaffold and gap extension** of haploid assemblies relative to a reference genome  
- **High-confidence haplotype read assignment** based on mapping quality, consistency, and error thresholds  
- **SNP filtering** with alignment context, base quality, read coverage, and repeat masking  
- **Identification of recombination events** (crossovers and gene conversions) through haplotype-consistent read support  
- **Event classification and alignment outputs** for downstream validation  

This makes the workflow broadly useful for both:  
- **Meiotic recombination studies** (e.g., sperm sequencing to map crossover landscapes)  
- **Somatic genome instability research** (e.g., recombination signatures in cancer genomes)  

---

## Workflow Overview

1. **Scaffold** haploid contig-level assemblies based on the reference assembly  
2. **Extend gaps (N’s)** to 30 kb (default) to avoid misinterpreting cross-gap mappings as crossovers  
3. **Align HiFi reads** to hap1 and hap2 assemblies separately  
4. **Filter primary reads** with:  
   - Unique mappings, MAPQ = 60  
   - Consistent alignment to the same chromosome and strand across haplotypes  
5. **Flag problematic reads** with:  
   - \>200 bp soft clipping (default), or  
   - \>100 sequence errors (mismatches to both haplotypes)  
6. **Filter SNPs** between haplotypes:  
   - At least 10 bp matched alignment flanking both sides  
   - \>400 bp away from read ends (Revio-specific)  
   - Base quality ≥ 40  
   - Outside low complexity or tandem repeat regions  
7. **Assign haplotype** to each retained HiFi read  
8. **Retain SNPs** only if covered by ≥3 reads with >95% haplotype consistency  
9. **Classify events** (crossover vs. gene conversion)  
10. **Align candidate reads** back to the reference assembly for reporting  

---

## Environment Requirements

- [Nextflow](https://github.com/nextflow-io/nextflow)  
- Singularity/Apptainer or Conda (see conda_environment.yml)
---

## Input Requirements

- Two **contig-level haploid assemblies** (hap1 and hap2)  
- **HiFi sequencing reads** (after QC)  
- A **chromosome-level reference assembly**

---

## Preparing Input Files

1. **Trim HiFi reads** using [HiFiAdapterFilt](https://github.com/sheinasim-USDA/HiFiAdapterFilt.git)  
2. **Assemble haplotypes** with [hifiasm](https://github.com/chhylp123/hifiasm.git)  
3. (Optional) **Validate assemblies** with BUSCO  
4. **Remove duplicates** using [purge_dups](https://github.com/dfguan/purge_dups.git)  

---

## Tips for Generating Balanced Haploid Assemblies

- Aim for **similar BUSCO scores** between hap1 and hap2; avoid imbalance where possible.  
- Fine-tune **hifiasm parameters** (`-s` and `--hom-cov`) to reduce collapse.  
- Alternatively, use **purge_dups** iteratively:  
  - Deduplicate hap1 and move duplicates into hap2  
  - Re-deduplicate hap2  
- This approach can **rescue over-collapsed contigs** and balance assembly completeness.


## Running pipeline

```bash
nextflow run main.nf \
  --hifi input_hifi \
  --ref ref_assembly.fa \
  --hap1 haploid_assembly_1.fa \
  --hap2 haploid_assembly_2.fa \
  --threads <num_of_cpus> \
  --profile standard \
  -c nextflow.config
```
All the following files and folders **must be in the same working directory** for the pipeline to run correctly:
```bash
project_directory/
├── hifi_sperm_recomb.sif # Singularity image file
├── main.nf               # Nextflow main workflow script
├── modules/              # Workflow modules
├── nextflow.config       # Nextflow configuration file
└── scripts/              # Custom scripts
```


## Parameters

| Parameter         | Description                                                                 | Default / Required |
|-------------------|-----------------------------------------------------------------------------|--------------------|
| `--hifi`          | HiFi sequencing reads input file (FASTQ/FASTA after QC).                    | **Required**       |
| `--ref`           | Chromosome-level reference assembly.                                        | **Required**       |
| `--hap1`          | Haploid assembly 1 (hap1).                                                  | **Required**       |
| `--hap2`          | Haploid assembly 2 (hap2).                                                  | **Required**       |
| `--leng_N_Gap`    | Length (in bp) to extend gaps (N’s) in assemblies to avoid cross-gap mapping.| `30000`            |
| `--soft_clip_num` | Maximum allowed soft-clipped bases per read; reads exceeding this are discarded. | `200`          |
| `--threads`       | Number of CPU threads to use.                                               | `8`                |
| `--help`          | Show the help message and exit.                                             | —                  |

---

## Profiles

Execution profiles define how and where the pipeline runs:

- **`-profile standard`**  
  Run locally on the current machine. Suitable for small datasets or testing.

- **`-profile slurm`**  
  Run on an HPC cluster with **SLURM** workload manager. Suitable for large datasets or production runs.


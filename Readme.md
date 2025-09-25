# Identifying Cross Overs and Gene Covnerstions from gamete HIFI sequening (Revio platform).
This script is designed to identify crossovers and gene conversion events from HiFi sequencing data (PacBio Revio platform) by comparing phased haploid assemblies and read alignments. It implements a series of stringent filtering steps on read mapping and SNP calling to ensure accurate detection of recombination signals.

While originally developed for gamete sequencing (e.g. sperm), the workflow is generalizable and can also be applied to bulk sequencing data, such as cancer tumor samples, where distinguishing haplotype-specific variants and detecting recombination-like events is important for studying genomic instability and clonal evolution.

Key features include:

Scaffold and gap-extension of haploid assemblies relative to a reference genome.

High-confidence read assignment to haplotypes, based on mapping quality, consistency, and error thresholds.

SNP filtering based on alignment context, base quality, read coverage, and repeat masking.

Identification of crossovers and gene conversions through haplotype-consistent read support.

Generation of event classifications and alignments for downstream validation.

This makes the script useful both for meiotic recombination studies (e.g. sperm sequencing to map crossover landscapes) and for somatic genome instability research (e.g. detecting recombination signatures in cancer genomes).


## What is the script doing?
    
    1. Scaffolding contig-level haploid assemblies based on the reference assembly.
    
    2. Extending N gaps to 30k (by default) to prevent cross-gap mappings that could be misinterpreted as crossovers.
    
    3. Align HIFI reads to haploid assemblies separately.
    
    4. Keep primary reads with a MAPQ = 60 that mapped uniquely and consistently to both haplotypes on the same chromosome and strand.
    
    5. Keep reads with >200 bp (by default) soft clipping or >100 sequence errors (sites mismatching both haplotypes).
    
    6. Further keep SNP (betwee haplotypes) with following criteria:
        6.1 Surrounded by at least 10bp of matched alignment on both sides of the SNP
        6.2 outside 400bp of either read ends for the Revio data
        6.3 Has a baseQ 40
        6.4 not located in low complexity and tandem repeats

    7. Assigning haplotype to each retained HIFI reads.

    8. Further retaining SNP to be covered by ≥3 reads with >95% consistency to the assigned haplotype.

    9. Classifying events.

    10. Aligning candidate reads to the reference assembly and output the result.





## Input file requirments:
    1. Two contig-level haploid assemblies (hap1 and hap2)
    2. HIFI sequencing reads after QC
    3. A reference assembly (chromosomal-level)


## Steps to generate input files:

    1. HIFI reads trimming with HiFiAdapterFilt (from GitHub https://github.com/sheinasim-USDA/HiFiAdapterFilt.git).

    2. HIFI assemblies generating using hifiasm (from GitHub https://github.com/chhylp123/hifiasm.git).
    When necessary:
    3. Haploid assemblies validation with BUSCO

    4. Duplication purging using purge_dups (GitHub https://github.com/dfguan/purge_dups.git)

Tipps:
Make sure the BUSCO scores of the two haploid assemblies are close to each other. Try to avoid imbalance as much as possible.
You can fine-tuning the hifiasm parameter (-s and --hom-cov) or using purge_dups to dedup hap1 and move the duplicates into hap2, then re-dedup hap2.
This will rescue some over-collapsed contigs and balance the completness.




    

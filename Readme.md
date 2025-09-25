# Identifying Cross Overs and Gene Covnerstions from gamete HIFI sequening (Revio platform).


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




    

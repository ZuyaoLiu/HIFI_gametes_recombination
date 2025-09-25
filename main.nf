#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

//Importing modules
include { PREPARING_HAPS } from './modules/preparing_haps.nf'
include { ALIGNING_AND_FILTERING } from './modules/aligning_and_filtering.nf'
include { ASSIGNING_HAPLOTYPES } from './modules/assigning_haplotypes.nf'
include { FILTERING_SNP } from './modules/filtering_SNP.nf'
include { IDENTIFYING_CANDIDATES } from './modules/identifying_candidates.nf'
include { MAPPING_TO_REF } from './modules/mapping_to_ref.nf'

params.help = false
if (params.help) {
    log.info """
================================================================================
 HiFi-Recomb pipeline
 Author : Zuyao Liu
 Date: 2025/09/04
 Version: 0.1.0

 Usage:
    nextflow run main.nf --samples samples.csv --ref ref.fa --outdir results

 Parameters:
    --hifi              HIFI reads input (required)
    --ref               Chromosome-level reference assembly (required)
    --hap1              Hap1 assembly (required)
    --hap2              Hap2 assembly (required)
    --leng_N_Gap        Length to extend (default: 30k)
    --soft_clip_num     Maximum soft-clipped bases per read (default: 200); reads exceeding this threshold are discarded.
    --threads           Number of threads to use (default: 8)  
    --help              Show this help page


 Profiles:
   -profile standard   Run locally
   -profile slurm      Run on HPC cluster
   -profile docker     Run with Docker/Singularity
================================================================================
"""
    exit 0
}


// Parameters chekcing
if( !params.hifi ) exit 1, "Missing required parameter: --hifi"
if( !params.ref  ) exit 1, "Missing required parameter: --ref"
if( !params.hap1 ) exit 1, "Missing required parameter: --hap1"
if( !params.hap2 ) exit 1, "Missing required parameter: --hap2"

if( !file(params.hifi).exists() ) exit 1, "HIFI file not found: ${params.hifi}"
if( !file(params.ref).exists()  ) exit 1, "Reference file not found: ${params.ref}"
if( !file(params.hap1).exists() ) exit 1, "Hap1 file not found: ${params.hap1}"
if( !file(params.hap2).exists() ) exit 1, "Hap2 file not found: ${params.hap2}"


workflow{

   // Setting input channels
   hifi_ch              = Channel.value( file(params.hifi) )
   ref_ch               = Channel.value( file(params.ref) )
   
   haps_ch = Channel.of(
      tuple('hap1', file( params.hap1 )),
      tuple('hap2', file( params.hap2 ))
      )

   leng_N_Gap_ch   = Channel.value( params.leng_N_Gap )
   soft_clip_num_ch = Channel.value( params.soft_clip_num )


   //Seting script channels
   extendN_py_ch = Channel.value( file('scripts/extending_N_Gaps.py') )

   filter_py_ch = Channel.value( file('scripts/filter_reads.py') )

   assigning_haplotypes_py_ch = Channel.value( file('scripts/comparing_haplotype_alignments.py') )

   trf_bed_py_ch = Channel.value( file('scripts/dat2bed.py') )

   filter_identify_candidates_py_ch = Channel.value( file('scripts/filter_identify_candidates.py') )
   
   classify_events_py_ch = Channel.value( file('scripts/classify_events.py') )

   mapping_to_ref_py_ch = Channel.value( file('scripts/map_read_snps_to_ref.py') )   




   //Main process
   prep = PREPARING_HAPS(haps_ch, ref_ch, extendN_py_ch, leng_N_Gap_ch)
   
   aln_filt = ALIGNING_AND_FILTERING(prep, hifi_ch, filter_py_ch)
               
   reads_assign =  ASSIGNING_HAPLOTYPES(aln_filt, assigning_haplotypes_py_ch, soft_clip_num_ch)
   
   snp_filtered = FILTERING_SNP(prep, reads_assign, trf_bed_py_ch)

   events = IDENTIFYING_CANDIDATES(snp_filtered, filter_identify_candidates_py_ch, classify_events_py_ch)

   map_to_ref = MAPPING_TO_REF(events, hifi_ch, ref_ch, mapping_to_ref_py_ch)
   
}


#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process EXTRACTING_CANDIDATES{
     label "mapping_to_ref"
    
    input:
        path(events)
        path(hifi_ch)
   
    output:
        path("candidates.fastq")

    script:

    """
    seqkit grep -f <(cat ${events} |grep -v "read_id"|cut -f1|sort -k1,1V|uniq)  ${hifi_ch} > candidates.fastq
    """

}

process MAPPING{
     label "mapping_to_ref"
    
    input:
        path(fastq)
        path(ref_ch)
   
    output:
        path("candidates_to_ref.bam")

    script:

    """
    set -euo pipefail
    minimap2 \
        -ax map-hifi --cs=short --eqx --MD \
        -t ${task.cpus} \
        ${ref_ch} \
        ${fastq}| \
    samtools sort -o candidates_to_ref.bam
    """
    
}

process IDENTIFYING_SITES{
    label "mapping_to_ref"
    publishDir "results/final_result", mode: 'copy'
    input:
        path(events)
        path(bam)
        path(mapping_to_ref_py_ch)
   
    output:
        path("Candidate_reads_SNP.onRef.txt")

    script:

    """
    python ${mapping_to_ref_py_ch} --bam ${bam} --tsv ${events} --out Candidate_reads_SNP.onRef.txt
    """
}
    


workflow MAPPING_TO_REF{
    take:
        events
        hifi_ch
        ref_ch
        mapping_to_ref_py_ch
    main:
        candidates_reads = EXTRACTING_CANDIDATES(events, hifi_ch)
        bam = MAPPING(candidates_reads, ref_ch)
        result = IDENTIFYING_SITES(events, bam, mapping_to_ref_py_ch)

    emit:
        result
}
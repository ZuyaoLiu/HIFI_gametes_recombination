#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process ALIGNING{
    label "aligning_and_filtering"
    tag { hap_name }

    input:
        tuple val(hap_name), path(hap)
        path(HIFI)
    
    output:
        tuple val(hap_name), path("${hap_name}.mq60.bam"), path("${hap_name}.mq60.bam.bai"), emit:bam

    script:

    """
    set -euo pipefail
    minimap2 \
        -ax map-hifi --cs=short --eqx --MD \
        -t ${task.cpus} \
        ${hap} \
        ${HIFI}| \
    samtools view -h -F 0x900 -q 60 -b |\
    samtools sort -o ${hap_name}.mq60.bam


  samtools index ${hap_name}.mq60.bam
    """
}



process FILTERING{
    label "aligning_and_filtering"

    publishDir "results/filtering_reads", mode: 'copy'

    input:
        tuple path(bam1), path(bai_1), path(bam2), path(bai_2)
        path(filter_py_ch)
    
    output:
        tuple path("hap1.mq60.filt.bam"), path("hap1.mq60.filt.bam.bai"), path("hap2.mq60.filt.bam"), path("hap2.mq60.filt.bam.bai"), path("retained_reads.list"), emit:bam

    script:

    """
    set -euo pipefail
    python ${filter_py_ch} \
        ${bam1} ${bam2} \
        -o retained_reads.list

    for i in hap1 hap2
    do
        samtools view -N retained_reads.list -b \${i}.mq60.bam | \
        samtools sort -o \${i}.mq60.filt.bam
        
        samtools index \${i}.mq60.filt.bam
    done
    """
}

workflow ALIGNING_AND_FILTERING{
    
    take:
        input_ch
        hifi_ch
        filter_py_ch
        
    
    main:
        aln_ch     = ALIGNING(input_ch, hifi_ch)
                        .collect(flat: false)
                        .map{ tuple1, tuple2 -> tuple tuple1[1], tuple1[2], tuple2[1], tuple2[2]  }
                          
        filt_ch    = FILTERING(aln_ch, filter_py_ch)

    emit:
        filt_ch

}



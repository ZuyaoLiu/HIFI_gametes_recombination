#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process SPLITTING_READS_LIST{
    label "assigning_haplotypes"
    
    input:
        path(retained_reads)

    
    output:
        path("*.reads.txt")

    script:

    """
    set -euo pipefail

    split -l 10000 -d -a 4 ${retained_reads} batch_

    for i in batch_*
    do
        mv "\$i" "\${i}.reads.txt"
    done

    """
}


process ASSIGNING_HAPS{
    label "assigning_haplotypes"
    

    input:
        
        tuple path(batch), path(bam1), path(bai_1), path(bam2), path(bai_2)
        path(assigning_haplotypes_py_ch)
        val(soft_clip_num)

    
    output:
        path ("${batch}.SNP.stats")

    script:
    def half_cpus = task.cpus.intdiv(2) 

    """
    set -euo pipefail

    samtools view -@ ${half_cpus} -N ${batch} ${bam1} -o ${batch}.hap1.bam &
    samtools view -@ ${half_cpus} -N ${batch} ${bam2} -o ${batch}.hap2.bam &
    wait

    python ${assigning_haplotypes_py_ch} \
        ${batch}.hap1.bam \
        ${batch}.hap2.bam \
        1> ${batch}.SNP.stats

    rm ${batch}.hap1.bam ${batch}.hap2.bam    
    """
}


process MERGE_READS_STATS{
    label "assigning_haplotypes"
    publishDir "results/reads_haplotype_assignment", mode: 'copy'
    input:
        path(batch_file)

    
    output:
        path("SNP.stats")

    script:

    """
    set -euo pipefail

    cat ${batch_file.join(' ')}  > temp
    head temp -n 1 > header
    grep -v "read_id" temp > stats
    cat header stats > SNP.stats
    """
}


workflow ASSIGNING_HAPLOTYPES{
    
    take:
        input_ch
        assigning_haplotypes_py_ch
        soft_clip_num_ch

        
    
    main:
        retained_reads = input_ch.map { it[4] } 
        split          = SPLITTING_READS_LIST(retained_reads)
        input_ch_mutant= input_ch.map { bam1, bai_1, bam2, bai_2, retained_reads -> tuple (bam1, bai_1, bam2, bai_2)} 
        assigned       = ASSIGNING_HAPS(split.flatten().combine(input_ch_mutant), assigning_haplotypes_py_ch, soft_clip_num_ch) 
                        .collect()

        merged_stats   = MERGE_READS_STATS( assigned  )
    emit:
        merged_stats


}


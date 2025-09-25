#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Scaffolding using RagTag
process SCAFFOLDING{
    label "preparing_haps"
    tag { hap_name }

    input:
        tuple val(hap_name), path(hap)
        path(ref)
    
    output:
        tuple val(hap_name), path("output/ragtag.scaffold.fasta"), emit:fasta

    script:

    """
    set -euo pipefail

    ragtag.py scaffold \
        ${ref} \
        ${hap} \
        -t ${task.cpus} \
        -u -w \
        --aligner minimap2 \
        -o output 
    """
}


// Extend Gaps to prevent cross-gap mappings
process EXTEND_GAPS{
    label "preparing_haps"
    tag { hap_name  }

    input:
        tuple val(hap_name), path(scaf)
        path(extendN_py_ch)
        val(leng_N_Gap_ch)
    
    output:
        tuple val(hap_name), path("ragtag.scaffold.fasta.N*.fasta"), emit:fasta

    script:

    """
    set -euo pipefail

    python ${extendN_py_ch} \
        ${scaf} \
        --gapsize ${leng_N_Gap_ch}

    """
}


workflow PREPARING_HAPS{
    
    take:
        haps_ch
        ref_ch
        extendN_py_ch
        leng_N_Gap_ch
    
    main:
        scaf_ch = SCAFFOLDING(haps_ch, ref_ch)
        ext_ch  = EXTEND_GAPS(scaf_ch, extendN_py_ch, leng_N_Gap_ch)

    emit:
        ext_ch

}

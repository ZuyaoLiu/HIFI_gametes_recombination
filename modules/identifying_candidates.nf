#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process SPLIT_CHR{
    label "identify_candidates"
    
    input:
        path(snp_filtered)
   
    output:
        path("*.SNP.filt.txt")

    script:

    """
    mkdir tmp
    head -n 1 ${snp_filtered} > header
    for i in \$(cat ${snp_filtered} |cut -f 5|grep -v "hap1_chr"|sort -k1,1V -T ./tmp |uniq)
    do
    awk -v x=\${i}  '\$5 == x' ${snp_filtered} > temp
    cat header temp > \${i}.SNP.filt.txt
    done
    """
}

process SNP_FILTER_FURTHER{
    label "identify_candidates"
    
    input:
        path(chr)
        path(filter_identify_candidates_py_ch)
    
    output:
        tuple path("${chr}_dropped.tsv"), path("${chr}_final_filt.tsv"), path("${chr}_trans_pairs.tsv") 

    script:
        
    """
    python ${filter_identify_candidates_py_ch} ${chr} ${chr}
    """
}

process MERGE_SNP_FILTER{
    label "identify_candidates"
    publishDir "results/identify_candidates", mode: 'copy'
    input:
        tuple path(dropped),
              path(filt),
              path(trans)
   
    output:
        tuple path("SNP.filt_dropped.tsv"), path("SNP.filt_final.tsv"), path("SNP.filt_trans_pairs.tsv")

    script:
     
    """
    cat ${dropped.join(' ')} > temp_dropped
    head -n 1 temp_dropped > header_dropped
    grep -v "read_id" temp_dropped > content_dropped
    cat header_dropped content_dropped > SNP.filt_dropped.tsv

    cat ${filt.join(' ')} > temp_filt
    head -n 1 temp_filt > header_filt
    grep -v "read_id" temp_filt > content_filt
    cat header_filt content_filt > SNP.filt_final.tsv

    cat ${trans.join(' ')} > temp_trans
    head -n 1 temp_trans > header_trans
    grep -v "read_id" temp_trans > content_trans
    cat header_trans content_trans > SNP.filt_trans_pairs.tsv


    """
}

process CLASSIFY_EVENTS{
    label "identify_candidates"
    
    input:
        path(final_filt)
        path(classify_events_py_ch)
   
    output:
        path("rec_events.tsv")

    script:
     
    """
    python ${classify_events_py_ch} ${final_filt}
    """
}

workflow IDENTIFYING_CANDIDATES{
    
    take:
        snp_filtered
        filter_identify_candidates_py_ch
        classify_events_py_ch

        
    
    main:
        chr           = SPLIT_CHR(snp_filtered).flatten() 
        snp_filtered_further        = SNP_FILTER_FURTHER( chr, filter_identify_candidates_py_ch)
                                        .collect(flat: false)
                                        .map { rows ->
                                        def dropped = rows.collect { it[0] }
                                        def filt    = rows.collect { it[1] }
                                        def trans   = rows.collect { it[2] }
                                        tuple(dropped, filt, trans)
                                        }
        merged_filt    = MERGE_SNP_FILTER(snp_filtered_further)
                       .map{dropped, final_filt, trans -> final_filt}
        
        events         = CLASSIFY_EVENTS(merged_filt, classify_events_py_ch) 
        
    emit:
        events 

}

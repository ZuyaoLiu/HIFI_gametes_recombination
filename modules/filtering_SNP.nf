#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process SDUST{
     label "filtering_SNP"
     tag { hap_name }
    input:
        tuple val(hap_name), path(fasta)
        
    
    output:
        tuple val(hap_name), path("${hap_name}.sdust")
        

    script:

    """
    sdust ${fasta} > ${hap_name}.sdust

    """
}

process TRF{
     label "filtering_SNP"
     tag { hap_name }
    input:
        tuple val(hap_name), path(fasta), path(script)
    
    output:
        tuple val(hap_name), path("${hap_name}.trf.bed")
       

    script:

    """
    #trf ${fasta} 2 6 6 80 10 50 500 -ngs -h -l 10 > ${hap_name}.trf
    #python ${script} --dat ${hap_name}.trf --out ${hap_name}.trf.bed --tool lobstr
    touch ${hap_name}.trf.bed
    """
}

process MERGE_BED{
     label "filtering_SNP"
     
    input:
        tuple val(hap_name), path(sdust), path(trf)
    
    output:
        tuple val(hap_name), path("${hap_name}.repeats.bed")

    script:

    """
    bedtools intersect -a ${sdust} -b ${trf} > ${hap_name}.repeats.bed
    """
}

process SNP_FILTER{
    label "filtering_SNP"
    publishDir "results/filtering_SNP", mode: 'copy'
    
    input:
        tuple path(hap1_repeats), path(hap2_repeats)
        path(reads_assign)
    
    output:
        path("SNP.filt.txt")

    script:

    """
    #Keeping SNP:
    #1. surrounded by at least 10bp of matched alignment on both sides of the SNP
    #2. SNPs that were outsied 400bp of either read ends for the Revio data
    #3. has a baseQ 40
    #4. not located in sdust and trf beds

    cat ${reads_assign} | head -n 1 > header
    cat ${reads_assign} | awk ' \$18 == "True" && \$4 == 40 && \$3 >=400 && (\$2-\$3+1) >= 400 ' > SNP.step123.filt

    cat SNP.step123.filt |cut -f 5,6|awk '{print \$1"\t"\$2-1"\t"\$2}'|sort -k1,1V -k2,2V -u > hap1.list
    cat SNP.step123.filt |cut -f 8,9|awk '{print \$1"\t"\$2-1"\t"\$2}'|sort -k1,1V -k2,2V -u > hap2.list

    bedtools intersect -a hap1.list -b ${hap1_repeats} -v |cut  -f 1,3 > hap1.snp.kept
    bedtools intersect -a hap2.list -b ${hap2_repeats} -v |cut  -f 1,3 > hap2.snp.kept

    awk 'NR==FNR {k[\$1 FS \$2]; next} (\$5 FS \$6) in k' hap1.snp.kept SNP.step123.filt > SNP.step1234.filt.hap1
    awk 'NR==FNR {k[\$1 FS \$2]; next} (\$8 FS \$9) in k' hap2.snp.kept SNP.step123.filt > SNP.step1234.filt.hap2

    mkdir tmp
    cat SNP.step1234.filt.hap1 SNP.step1234.filt.hap2 | \
        sort -k1,1V -k3,3V \
        -S 8G \
        --parallel=${task.cpus} \
        -T ./tmp | \
        uniq -d | awk ' \$5 == \$8 ' > 1   

    cat header 1  > SNP.filt.txt
    """
}

workflow FILTERING_SNP{
    
    take:
        prep
        reads_assign
        trf_bed_py_ch

    main:
        sdust_input = prep
        trf_input =  prep.combine(trf_bed_py_ch)
        sdust_mask = SDUST(sdust_input)
        trf_mask = TRF( trf_input )
        repeat_merge_input_ch = sdust_mask.join(trf_mask)
        repeat_merge = MERGE_BED(repeat_merge_input_ch)
                        .collect()
                        .map{ hap1, bed1, hap2, bed2 -> tuple bed1, bed2}
        snp_filter = SNP_FILTER(repeat_merge, reads_assign)
       
    emit:
        snp_filter
}
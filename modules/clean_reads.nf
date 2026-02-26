process clean_reads {
    // Clean reads (adapter and read length filter) with fastp
    
    tag '1'
    publishDir params.outDir + '/QC', mode: 'copy'
    
    input:
      tuple val(params.sampleid), path(reads)

    output:
      path("${reads[0].simpleName}_fastp.fastq.gz"), emit: read1
      path("${reads[1].simpleName}_fastp.fastq.gz"), emit: read2
      path("${params.sampleid}.fastp.html")
      path("${params.sampleid}_fastpmerged.fastq.gz"), emit: merged
      path("${params.sampleid}.fastp.json")

    script:

      // def (read1, read2) = reads

      """
        fastp -i ${reads[0]} -I ${reads[1]} \
        -o ${reads[0].simpleName}_fastp.fastq.gz -O ${reads[1].simpleName}_fastp.fastq.gz \
        --merge --merged_out ${params.sampleid}_fastpmerged.fastq.gz \
        --overlap_diff_limit 0 \
        --trim_poly_x --trim_poly_g --length_required 100 --thread ${params.threads} --detect_adapter_for_pe\
        --json ${params.sampleid}.fastp.json --html ${params.sampleid}.fastp.html
      """
}
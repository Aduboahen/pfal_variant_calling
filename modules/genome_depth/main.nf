process genome_depth {
    // Genome depth with mosdepth 

    tag '4'
    publishDir params.outDir + '/QC', mode: 'copy'

    input:
      path(bam)
      path(bam_index)

    output:
      path "${params.sampleid}*"

    script:
    """
      mosdepth --threads ${params.threads} ${params.sampleid} ${bam}
    """
}

process genome_depth_ampseq {
    // Genome depth with mosdepth 

    tag '4'
    publishDir params.outDir + '/QC', mode: 'copy'

    input:
      path(bam)
      path(bam_index)

    output:
      path "${params.sampleid}*"

    script:
    """
      mosdepth --threads ${params.threads} ${params.sampleid} ${bam}
      samtools bedcov ${params.bedfile} ${bam} > ${params.sampleid}.basecount.tsv
    """
}
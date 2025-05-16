process variant_calling{
    // call variants with bcftools

    tag '5'
    publishDir params.outDir + "/variants", mode: 'copy'

    input:
      path(bam)
      path(bam_index)

    output:
      path("${params.sampleid}.vcf.gz")

    script:
    """
      bcftools mpileup -a DP -B -O u -m 4 --threads ${params.threads} -f \
      ${params.reference} ${bam}\
      | bcftools call -mv -O z -o ${params.sampleid}.vcf.gz
    """
}


// call variants with bcftools for ampseq

process variant_calling_ampseq{
    tag '5'

    publishDir params.outDir + "/variants", mode: 'copy'

    input:
      path(bam)
      path(bam_index)

    output:
      path("${params.sampleid}.vcf.gz")

    script:
    """
      bcftools mpileup -a DP -B -O u -m 4 --threads ${params.threads} -f \
      -R ${params.bed_file} \
      ${params.reference} ${bam}\
      | bcftools call -mv -O z -o ${params.sampleid}.vcf.gz
    """
}
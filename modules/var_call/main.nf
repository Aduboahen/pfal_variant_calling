process variant_calling{
    // call variants with bcftools

    tag '5'
    publishDir params.outDir + "/variants", mode: 'copy'

    input:
      path(bam)
      path(bam_index)

    output:
      path("${params.sampleid}.vcf")

    script:
    """
      bcftools mpileup -a 'FORMAT/AD,FORMAT/DP' -B -O v -m 4 --threads ${params.threads} -f \
      ${params.reference} ${bam} \
      | bcftools call -mv -O z -o ${params.sampleid}.vcf
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
      path("${params.sampleid}.bcftools.vcf")

    script:
    """
      bcftools mpileup -a 'FORMAT/AD,FORMAT/DP' --skip-indels \
      --threads ${params.threads} -f ${params.reference} ${bam} \
      -R ${params.bedfile}\
      | bcftools call -m -V indels \
      -Ov -o ${params.sampleid}.bcftools.vcf
    """
}



// call variants with lofreq for ampseq

  process variant_calling_lofreq{
    tag '5'

    publishDir params.outDir + "/variants", mode: 'copy'

    input:
      path(bam)
      path(bam_index)

    output:
      path("${params.sampleid}.lofreq.vcf")

    script:
    """
      lofreq call-parallel --pp-threads ${params.threads} \
      -f ${params.reference} -l ${params.bedfile} \
      -o ${params.sampleid}.lofreq.vcf ${bam}
    """
}


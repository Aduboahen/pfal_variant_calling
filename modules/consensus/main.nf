process consensus{
  // This process generates a consensus sequence from the VCF file and the reference genome.
  // It uses bcftools to create a consensus FASTA file based on the VCF and reference genome.
  // The output is a FASTA file containing the consensus sequence.
  // The process also modifies the header of the FASTA file to include the sample ID.
  // The input includes the VCF file and a bed file for non-covered regions.
  // The output is a FASTA file with the consensus sequence.

    tag '10'
    publishDir params.outDir + "/consensus", mode: 'copy'

    input:
      path(vcf)
      path(noncov)

    output:
      path("${params.sampleid}.consensus.fasta")

    script:
        """
        bcftools view -O z -o ${params.sampleid}.vcf.gz ${vcf}
  
        bcftools index ${params.sampleid}.vcf.gz
  
        bcftools consensus -f ${params.reference}\
        --mark-snv lc\
        --mark-del -\
        --mask ${noncov}\
      ${params.sampleid}.vcf.gz  > ${params.sampleid}.consensus.fasta

      sed -i 's/>Pf3D7/>${params.sampleid}/;s/_v3/''/g' ${params.sampleid}.consensus.fasta
      """
}
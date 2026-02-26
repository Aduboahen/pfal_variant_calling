process ubiquitous_variant {

  // Process 3C: ubiquitous_variant                                                         
  // Find the variants where we are uncertain that the variation is correct
  // so that they can be replaced/masked in the consensus seq by N
    publishDir params.outDir + '/uncovered', mode: 'copy'
    input:
        path vcf
    output:
        path "${params.sampleid}_ubiq.bed"
    script:
        """
        bcftools query -f "%CHROM\\t%POS\\t%END\\n" -i \
        "(%QUAL<${params.varqual} || \
        %MAX(FORMAT/AD[0:1])/%MAX(FORMAT/DP)<${params.varthres}) && \
        INFO/DP>=${params.depth}" ${vcf} | awk '{{print(\$1 \"\\t\" \$2 \"\\t\" \$3 \    "\\tubiquitous_variant\")}}' > ${params.sampleid}_ubiq.bed
        """
} 


process non_covered_regions{
  // non-covered regions
  // using bedtools genomecov
  // and bedtools subtract
  // to get regions of low coverage
  // and remove them from the vcf file

    publishDir params.outDir + '/uncovered', mode: 'copy'

    input:
      path(bam)
      path(vcf)

    output:
      path("${params.sampleid}_noncov.bed")

    script:
    """
      bedtools genomecov -ibam ${bam} -bga | \
      awk '\$4 < ${params.depth}' | \
      awk '{{print(\$1 \"\\t\" \$2 + 1 \"\\t\" \$3 \"\\tlow_coverage\")}}' |\
      bedtools subtract -a - -b ${vcf} > ${params.sampleid}_noncov.bed
    """
}
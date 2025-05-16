process non_covered_regions{
  // non-covered regions
  // using bedtools genomecov
  // and bedtools subtract
  // to get regions of low coverage
  // and remove them from the vcf file

    tag '9'
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
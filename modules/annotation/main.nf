process annotation_snps {
    // Annotate variants with snpEff
    // 
    // This process annotates the variants in the VCF file using snpEff and generates a summary HTML file, a CSV file with statistics, and a table of annotations.
    // The input is a VCF file, and the output includes the annotated VCF file, summary HTML, CSV statistics, and annotation table.
    // The process uses the params.vcf2table script to convert the annotated VCF into a tabular format.
    //
    // Input:
    // - vcf: The input VCF file containing variants.
    //
    // Output:
    // - snpEffstats: A CSV file containing statistics about the annotation.
    // - summary.html: An HTML summary of the annotation.
    // - genes.txt: A text file containing gene information.
    // - annot.table.tsv: A tab-separated values (TSV) file containing the annotation table.
    // - snpEff.ann.vcf: The annotated VCF file.

    tag '8a'
    publishDir params.outDir + "/annotation", mode: 'copy'
    // publishDir params.outDir + '/QC', mode: 'copy', pattern: "${params.sampleid}.snpEff.snps.csv"

    input:
      path(vcf)

    output:
      path("${params.sampleid}.snpEff.snps.csv"), emit: snpEffstats
      path("${params.sampleid}.summary.snps.html")
      path("${params.sampleid}.snpEff.snps.genes.txt")
      path("${params.sampleid}.annot.table.snps.tsv")
      path("${params.sampleid}.snpEff.ann.snps.vcf")

    script:
    """
      snpEff ann -v -noShiftHgvs -ud 0 -strict -s -hgvs1LetterAa \
      ${params.sampleid}.summary.snps.html Plasmodium_falciparum\
      -csvStats ${params.sampleid}.snpEff.snps.csv ${vcf} > ${params.sampleid}.snpEff.ann.snps.vcf

      ${params.vcf2table} ${params.sampleid}.snpEff.ann.snps.vcf --sample ${params.sampleid} -ad -e -o ${params.sampleid}.annot.table.snps.tsv
    """
}


process annotation_indels {
  // This process is similar to the previous one but specifically for indels

    tag '8b'
    publishDir params.outDir + "/annotation", mode: 'copy'
    // publishDir params.outDir + '/QC', mode: 'copy', pattern: "${params.sampleid}.snpEff.indels.csv"

    input:
      path(vcf)

    output:
      path("${params.sampleid}.snpEff.indels.csv"), emit: snpEffstats
      path("${params.sampleid}.summary.indels.html")
      path("${params.sampleid}.snpEff.indels.genes.txt")
      path("${params.sampleid}.annot.table.indels.tsv")
      path("${params.sampleid}.snpEff.ann.indels.vcf")

    script:
    """
      snpEff ann -v -noShiftHgvs -ud 0 -strict -s -hgvs1LetterAa  \
      ${params.sampleid}.summary.indels.html Plasmodium_falciparum\
      -csvStats ${params.sampleid}.snpEff.indels.csv ${vcf} > ${params.sampleid}.snpEff.ann.indels.vcf

      ${params.vcf2table} ${params.sampleid}.snpEff.ann.indels.vcf --sample ${params.sampleid} -ad -e -o ${params.sampleid}.annot.table.indels.tsv
    """
}

process annotation {
    // Annotate variants with snpEff
    // 
    // This process annotates the variants in the VCF file using snpEff and generates a summary HTML file, a CSV file with statistics, and a table of annotations.
    // The input is a VCF file, and the output includes the annotated VCF file, summary HTML, CSV statistics, and annotation table.
    // The process uses the params.vcf2table script to convert the annotated VCF into a tabular format.
    //
    // Input:
    // - vcf: The input VCF file containing variants.
    //
    // Output:
    // - snpEffstats: A CSV file containing statistics about the annotation.
    // - summary.html: An HTML summary of the annotation.
    // - genes.txt: A text file containing gene information.
    // - annot.table.tsv: A tab-separated values (TSV) file containing the annotation table.
    // - snpEff.ann.vcf: The annotated VCF file.

    tag '8'
    publishDir params.outDir + "/annotation", mode: 'copy'
    publishDir params.outDir + '/QC', mode: 'copy', pattern: "${params.sampleid}.snpEff.csv"

    input:
      path(vcf)

    output:
      path("${params.sampleid}.snpEff.csv"), emit: snpEffstats
      path("${params.sampleid}.summary.html")
      path("${params.sampleid}.snpEff.genes.txt")
      path("${params.sampleid}.annot.table.tsv")
      path("${params.sampleid}.snpEff.ann.vcf")

    script:
    """
      snpEff ann -v -noShiftHgvs -ud 0 -strict -s -hgvs1LetterAa \
      ${params.sampleid}.summary.html Plasmodium_falciparum\
      -csvStats ${params.sampleid}.snpEff.csv ${vcf} > ${params.sampleid}.snpEff.ann.vcf

      ${params.vcf2table} ${params.sampleid}.snpEff.ann.vcf --sample ${params.sampleid} -ad -e -o ${params.sampleid}.annot.table.tsv
    """
}
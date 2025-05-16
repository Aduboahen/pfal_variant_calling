

process multiQC {
  // MultiQC report generation
  // This process generates a MultiQC report from the output of various tools
  // used in the workflow. It collects the results from different processes and
  // creates a comprehensive report in the specified output directory.
  // The report includes information from tools like BCFtools, SnpEff, and
  // Mosdepth, providing an overview of the quality control metrics and
  // statistics for the variant calling process.

  tag '10'
  publishDir params.outDir + '/QC', mode: 'copy'

  input:
    path(snpEffStats)

  output:
    path "*"
    path(".command.*")

  script:
        """
        multiqc ${params.outDir}/QC
        """
}
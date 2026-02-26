#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


include { mapping_ampseq } from '../modules/mapping'
include { genome_depth_ampseq } from '../modules/genome_depth'
include { variant_calling_ampseq ; variant_calling_lofreq } from '../modules/variant_calling'
include { filter_snps ; filter_lofreq } from '../modules/filtering'
include { annotation_bcftools ; annotation_lofreq } from '../modules/annotation'
// include { multiQC } from '../modules/multiQC'


process run_parameters {
  // This process generates a parameters file that includes the operational
  // parameters used in the workflow.

    publishDir params.outDir, mode: 'copy'

    output:
        path("parameters.txt")

    script:
        """
        touch parameters.tsv
        echo "Parameter\tValue" >> parameters.txt
        echo "Mapping quality:\t>=${params.mapqual}" >> parameters.txt
        echo "Mutation frequency:\t>=${params.varthres}" >> parameters.txt
        echo "Variant quality:\t>${params.varqual}" >> parameters.txt
        echo "Minimum sequence depth:\t${params.depth}" >> parameters.txt
        echo "Minimum strand bias score:\t${params.strandthres}" >> parameters.txt
        echo "BWA seed length:\t${params.seed_length}" >> parameters.txt
        """
}


workflow {

  // channel to get reads as tuples


  reads_ch = channel.fromFilePairs(params.reads)
    .ifEmpty { "no such files" }

  log.info(
    """
P. falciparum variant calling

========Sources===============
cwd       : ${projectDir}
sample    : ${params.sampleid}
reads     : ${params.reads}
outDir    : ${params.outDir}
threads   : ${params.threads}

=======Variant filters========
variant quality   : ${params.varqual}
variant threshold : ${params.varthres}
variant depth     : ${params.depth}
mapping quality   : ${params.mapqual}

=======Reference=======
reference : ${params.reference}
bedfile   : ${params.bedfile}

=======Author=======
James Osei-Mensa
oseimensa@kccr.de
"""
  )
  // run the pipeline

  reads_ch = channel.fromFilePairs(params.reads).ifEmpty { error("No such files") }
  mapping_ampseq(reads_ch)
  genome_depth_ampseq(mapping_ampseq.out.bam, mapping_ampseq.out.bam_index)
  variant_calling_ampseq(mapping_ampseq.out.bam, mapping_ampseq.out.bam_index)
  filter_snps(variant_calling_ampseq.out)
  annotation_bcftools(filter_snps.out)
  variant_calling_lofreq(mapping_ampseq.out.bam, mapping_ampseq.out.bam_index)
  filter_lofreq(variant_calling_lofreq.out)
  annotation_lofreq(filter_lofreq.out)
  run_parameters()
}

#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// operational parameters

params.outDir = "${params.outDir}"
params.reads = "${projectDir}/test_fq/test_{1,2}.fastq.gz"
// params.read1    = "$projectDir/test_fq/test_1.fastq.gz"
// params.read2    = "$projectDir/test_fq/test_2.fastq.gz"
params.sampleid = "${params.sampleid}"

params.threads = 6
params.reference = "${projectDir}/reference/Pf3D7.fasta"
params.bedfile = "${projectDir}/reference/genes.bed"

// filtering parameters

params.varqual = 15
params.varthres = 0.4
params.mapqual = 15
params.depth = 5
params.seed_length = 100


// tools

params.vcf2table = "${projectDir}/scripts/vcf2table.py"
params.parse_stats = "${projectDir}/scripts/parse_stats.py"

// include { clean_reads } from './modules/clean_reads'
include { mapping_ampseq } from './modules/mapping'
include { genome_depth_ampseq } from './modules/genome_depth'
include { variant_calling_ampseq ; variant_calling_lofreq } from './modules/var_call'
include { filter_ampseq } from './modules/filter'
include { filter_lofreq } from './modules/filter'
include { annotation; annotation_filt; annotation_lofreq; annotation_lofreq2 } from './modules/annotation'
// include { multiQC } from './modules/multiQC'
include { run_parameters } from './modules/run_parameters'


workflow {

  // channel to get reads as tuples


  reads_ch = Channel.fromFilePairs(params.reads)
    .ifEmpty { "no such files" }

  log.info(
    """
P. falciparum variant calling

========Sources===============
working directory  : ${projectDir}
sample    : ${params.sampleid}
reads     : ${params.reads}
outDir    : ${params.outDir}
threads   : ${params.threads}


=======Variant filters========
variant quality   : ${params.varqual}
variant threshold  : ${params.varthres}
mapping quality   : ${params.mapqual}
read depth at variant site   : ${params.depth}


=======Reference=======
reference : ${params.reference}
bedfile   : ${params.bedfile}

=======Author=======
James Osei-Mensa
oseimensa@kccr.de
"""
  )
  // run the pipeline

  reads_ch = Channel.fromFilePairs(params.reads).ifEmpty { error("No such files") }
  // clean_reads(reads_ch)
  // mapping(clean_reads.out.read1, clean_reads.out.read2, clean_reads.out.merged)
  mapping_ampseq(reads_ch)
  genome_depth_ampseq(mapping_ampseq.out.bam, mapping_ampseq.out.bam_index)
  variant_calling_ampseq(mapping_ampseq.out.bam, mapping_ampseq.out.bam_index)
  annotation(variant_calling_ampseq.out)
  filter_ampseq(variant_calling_ampseq.out)
  annotation_filt(filter_ampseq.out)
  variant_calling_lofreq(mapping_ampseq.out.bam, mapping_ampseq.out.bam_index)
  annotation_lofreq2(variant_calling_lofreq.out)
  filter_lofreq(variant_calling_lofreq.out)
  annotation_lofreq(filter_lofreq.out)
  run_parameters()
}

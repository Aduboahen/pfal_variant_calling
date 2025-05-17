#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// operational parameters

params.outDir           = "$params.outDir"
params.params.reads     = "$projectDir/test_fq/test_{1,2}.fastq.gz"
params.sampleid         = "$params.sampleid"

params.threads = 6
params.reference        = "$projectDir/reference/Pf3D7.fasta"
params.bed_file				 = "$projectDir/reference/amplicons.bed"

// filtering parameters

params.varqual   = 15
params.varthres  = 0.80
params.mapqual   = 15 
params.depth     = 5


// tools

params.vcf2table        = "$projectDir/scripts/vcf2table.py"
params.parse_stats      = "$projectDir/scripts/parse_stats.py"

include { clean_reads } from './modules/clean_reads'
include { mapping } from './modules/mapping'
include { mark_duplicates } from './modules/mark_dup'
include { genome_depth } from './modules/genome_depth'
include { variant_calling_ampseq } from './modules/var_call'
include { filter; filter_snps; filter_indels; merge_vcf } from './modules/filter'
include { non_covered_regions } from './modules/non_cov_reg'
include { annotation; annotation_snps; annotation_indels } from './modules/annotation'
include { consensus } from './modules/consensus'
include { genome_stats } from './modules/genome_stats'
include { multiQC } from './modules/multiQC'
include { run_parameters } from './modules/run_parameters'


workflow {

// channel to get reads as tuples


reads_ch = Channel
    .fromFilePairs(params.reads)
    .ifEmpty{"no such files"}

log.info"""
P. falciparum variant calling

========Sources===============
codeBase  : "$projectDir"
sample    : "$params.sampleid"
reads     : "$params.reads"
outDir    : "$params.outDir"
threads   : 6


=======Variant filters========
varqual   : $params.varqual
varthres  : $params.varthres
mapqual   : $params.mapqual
depth     : $params.depth


=======Reference Pf3D7=======
reference : "$projectDir/reference/Pf3D7.fasta"


=======Author=======
James Osei-Mensa
oseimensa@kccr.de
"""
// run the pipeline
  
  reads_ch = Channel.fromFilePairs( params.reads ).ifEmpty { error "No such files" }
  clean_reads(reads_ch)
  mapping(clean_reads.out.read1, clean_reads.out.read2, clean_reads.out.merged)
  mark_duplicates(mapping.out.bam)
  genome_depth(mark_duplicates.out.bam, mark_duplicates.out.bam_index)
  variant_calling_ampseq(mark_duplicates.out.bam, mark_duplicates.out.bam_index)
  filter(variant_calling_ampseq.out)
  // filter_snps(variant_calling.out)
  // filter_indels(variant_calling.out)
  // merge_vcf(filter_snps.out.snps, filter_indels.out.indels, filter_snps.out.snps_index, filter_indels.out.indels_index)
  // annotation_snps(filter_snps.out.snps)
  // annotation_indels(filter_indels.out.indels)
  annotation(filter.out)
  non_covered_regions(mark_duplicates.out.bam, variant_calling_ampseq.out)
  consensus(variant_calling_ampseq.out, non_covered_regions.out)
  genome_stats(consensus.out)
  // // multiQC(annotation_indels.out)
  run_parameters()
}
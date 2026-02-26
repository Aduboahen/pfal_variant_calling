#!/usr/bin/env nextflow

nextflow.enable.dsl=2


include { clean_reads } from '../modules/clean_reads'
include { mapping } from '../modules/mapping'
include { mark_duplicates } from '../modules/mark_duplicates'
include { genome_depth } from '../modules/genome_depth'
include { variant_calling; variant_calling_lofreq } from '../modules/variant_calling'
include { filter; filter_snps; filter_indels; merge_vcf; filter_lofreq } from '../modules/filtering'
include { non_covered_regions } from '../modules/non_covered_regions'
include { annotation_bcftools; annotation_snps; annotation_indels; annotation_lofreq } from '../modules/annotation'
include { consensus } from '../modules/consensus'
include { genome_stats } from '../modules/genome_stats'
include { multiQC } from '../modules/multiQC'



process run_parameters {
  // This process generates a parameters file that includes the operational
  // parameters used in the workflow.

    publishDir params.outDir, mode: 'copy'

    output:
        path("parameters.txt")

    script:
        """
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

reads_ch = channel
    .fromFilePairs(params.reads)
    .ifEmpty{"no such files"}

log.info"""
P. falciparum variant calling

========Sources===============
codeBase  : $projectDir
sample    : $params.sampleid
reads     : $params.reads
outDir    : $params.outDir
threads   : 6


=======Variant filters========
varqual   : $params.varqual
varthres  : $params.varthres
mapqual   : $params.mapqual
depth     : $params.depth


=======Reference Pf3D7=======
reference : $params.reference


=======Author=======
James Osei-Mensa
oseimensa@kccr.de
"""
// run the pipeline
  
  reads_ch = channel.fromFilePairs( params.reads ).ifEmpty { error "No such files" }
  clean_reads(reads_ch)
  mapping(clean_reads.out.read1, clean_reads.out.read2, clean_reads.out.merged)
  mark_duplicates(mapping.out.bam)
  genome_depth(mark_duplicates.out.bam, mark_duplicates.out.bam_index)
  // variant_calling(mark_duplicates.out.bam, mark_duplicates.out.bam_index)
  variant_calling_lofreq(mark_duplicates.out.bam, mark_duplicates.out.bam_index)
  // filter(variant_calling.out)
  filter_lofreq(variant_calling_lofreq.out)
  // filter_snps(variant_calling.out)
  // filter_indels(variant_calling.out)
  // merge_vcf(filter_snps.out.snps, filter_indels.out.indels, filter_snps.out.snps_index, filter_indels.out.indels_index)
  // annotation_snps(filter_snps.out.snps)
  // annotation_indels(filter_indels.out.indels)
  // annotation(filter.out)
  non_covered_regions(mark_duplicates.out.bam, variant_calling_lofreq.out)
  annotation_lofreq(filter_lofreq.out)
  // consensus(variant_calling.out, non_covered_regions.out)
  consensus(variant_calling_lofreq.out, non_covered_regions.out)
  genome_stats(consensus.out)
  multiQC(annotation_bcftools.out.snpEffstats)
  run_parameters()
}
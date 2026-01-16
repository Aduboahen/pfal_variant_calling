#!/bin/bash
# coverage_analysis.sh

# Purpose: Script to run a Nextflow workflow to analyze coverage of amplicon sequencing data

# Synopsis: ./coverage_analysis.sh

# Positional arguments (defaults):
#   None

# Behavior:
#   - Iterates over files matching ./input/*R1*.fastq.gz
#   - Derives sample ID by stripping the R1 suffix and extension from the filename
#   - Invokes: nextflow run ./coverage.nf --sampleid <sampleid> --reads ./input/<sampleid>_{R1,R2}.fastq.gz --outDir ./output/coverage

# Requirements:
#   - nextflow available in PATH
#   - Input reads named consistently as <sampleid>_R1.fastq.gz and <sampleid>_R2.fastq.gz (or as configured)
#   - Reference files at ./reference/amplicons.bed and ./reference/Pf3D7.fasta (or adjust script)

# Notes:

workflow="${1:-./workflows/coverage.nf}"
input="${2:-./input/}"
suffix1="${3:-_R1}"
suffix2="${4:-.fastq.gz}"
thread="${5:-8}"

for file in ./input/*${suffix1}${suffix2}*; do
	sampleid=$(basename $file ${suffix1}${suffix2})

	nextflow run ${workflow} --sampleid ${sampleid} --reads ./input/"${sampleid}_{R1,R2}.fastq.gz" --outDir "./output/"coverage_${today}" --threads ${thread} -resume
done
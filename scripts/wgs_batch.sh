#!/bin/bash

today=$(date +"%d-%B-%y")
workflow="${1:-./BCFTOOLS_dsl2.nf}"
input="${2:-./input/}"
suffix1="${3:-_R1}"
suffix2="${4:-.fastq.gz}"
thread="${5:-6}"

for file in ./input/*${suffix1}${suffix2}*; do
	sampleid=$(basename $file ${suffix1}${suffix2})

	nextflow run ${workflow} --sampleid ${sampleid} --reads ./input/"${sampleid}_{R1,R2}.fastq.gz" --outDir ./output/${today} --threads ${thread}  --reference ./reference/Pf3D7.fasta -resume

done
# ampseq_batch.sh
# Purpose:
#   Batch-launch a Nextflow amplicon sequencing workflow for paired-end samples found in an input directory.
#
# Synopsis:
#   ./ampseq_batch.sh [WORKFLOW] [INPUT_DIR] [SUFFIX_R1] [SUFFIX_R2] [THREADS]
#
# Positional arguments (defaults):
#   WORKFLOW   - Path to Nextflow workflow script (default: ./ampseq.nf)
#   INPUT_DIR  - Directory containing input reads (default: ./input/)
#   SUFFIX_R1  - R1 filename suffix to match (default: _R1)
#   SUFFIX_R2  - R2 filename suffix/extension to match (default: .fastq.gz)
#   THREADS    - Number of threads to pass to workflow (default: 6)
#
# Behavior:
#   - Iterates over files matching INPUT_DIR/*${SUFFIX_R1}${SUFFIX_R2}*
#   - Derives sample ID by stripping the R1/R2 suffix and extension from the filename
#   - Invokes: nextflow run WORKFLOW --sampleid <sampleid> --reads ./input/"${sampleid}_{R1,R2}.fastq.gz"
#             --outDir ./output/<DATE> --threads THREADS --bedfile ./reference/amplicons.bed
#             --reference ./reference/Pf3D7.fasta -resume
#   - Uses a dated output directory (DAY-Month-YY)
#
# Requirements:
#   - nextflow available in PATH
#   - Input reads named consistently as <sampleid>_R1.fastq.gz and <sampleid>_R2.fastq.gz (or as configured)
#   - Reference files at ./reference/amplicons.bed and ./reference/Pf3D7.fasta (or adjust script)
#
# Notes:
#   - Adjust WORKFLOW and file paths as needed for your environment.
#   - The script appends -resume to retain Nextflow cache between runs.
#   - Ensure execute permission: chmod +x ampseq_batch.sh
#
# Example:
#   ./ampseq_batch.sh ./ampseq.nf ./input/ _R1 .fastq.gz 8
#!/bin/bash

today=$(date +"%d-%B-%y")
workflow="${1:-./workflows/ampseq.nf}"
input="${2:-./input/}"
suffix1="${3:-_R1}"
suffix2="${4:-.fastq.gz}"
thread="${5:-8}"

for file in ./input/*${suffix1}${suffix2}*; do
	sampleid=$(basename $file ${suffix1}${suffix2})

	nextflow run ${workflow} --sampleid ${sampleid} --reads ./input/"${sampleid}_{R1,R2}.fastq.gz" --outDir ./output/${today} --threads ${thread} -resume
done
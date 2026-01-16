# Batch runner for Nextflow whole-genome sequencing workflow
#
# Description:
#  Iterates over paired-end FASTQ files in the input directory and launches a Nextflow
#  run for each sample. Output for each run is written to ./output/<date>.
#
# Usage:
#  ./wgs_batch.sh [WORKFLOW] [INPUT_DIR] [SUFFIX_R1] [SUFFIX_R2] [THREADS]
#
# Positional arguments (defaults shown):
#  WORKFLOW   Path to Nextflow workflow script (default: ./workflows/wgs.nf)
#  INPUT_DIR  Directory containing input FASTQ files (default: ./input/)
#  SUFFIX_R1  R1 filename suffix to match (default: _R1)
#  SUFFIX_R2  R2 filename suffix/extension to match (default: .fastq.gz)
#  THREADS    Number of threads to pass to the workflow (default: 6)
#
# Expected input file naming:
#  Paired reads must be named as: <sample>_R1.fastq.gz and <sample>_R2.fastq.gz
#  The script searches for files matching: ./input/*${SUFFIX_R1}${SUFFIX_R2}*
#
# Behavior:
#  - For each matching R1 file, sample ID is derived by stripping SUFFIX_R1 and SUFFIX_R2.
#  - Invokes Nextflow:
#      nextflow run <WORKFLOW> --sampleid <sample> --reads ./input/"<sample>_{R1,R2}.fastq.gz"
#        --outDir ./output/<DD-Month-YY> --threads <THREADS> -resume
#  - Uses today's date (format: DD-Month-YY) to create the output subdirectory.
#
# Requirements & notes:
#  - Nextflow must be installed and available in PATH.
#  - Run the script from the repository root (or adjust relative paths).
#  - Filenames containing whitespace or unusual characters may cause issues; use safe names.
#  - The script uses -resume to continue previous Nextflow runs when possible.
#   - Reference files at ./reference/Pf3D7.fasta (or adjust script)

#
# Example:
#  ./wgs_batch.sh ./workflows/wgs.nf ./input/ _R1 .fastq.gz 8
#
# Author: Generated documentation
# Date: auto-generated
#!/bin/bash

today=$(date +"%d-%B-%y")
workflow="${1:-./workflows/wgs.nf}"
input="${2:-./input/}"
suffix1="${3:-_R1}"
suffix2="${4:-.fastq.gz}"
thread="${5:-6}"

for file in ./input/*${suffix1}${suffix2}*; do
	sampleid=$(basename $file ${suffix1}${suffix2})

	nextflow run ${workflow} --sampleid ${sampleid} --reads ./input/"${sampleid}_{R1,R2}.fastq.gz" --outDir ./output/${today} --threads ${thread} -resume

done
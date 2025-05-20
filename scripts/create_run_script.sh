#!/bin/bash

# bash script to generate a script (run_script_bs.sh) for running in batch mode
# the generated script contains one line for each sample
# the generated script must be run from within the docker container
# arguments:
# suffix: suffix of the input files (e.g. _R1)
# workflow: name of the workflow to run (e.g. BCFTOOLS.nf)
# threads: number of threads to use
# usage: ./create_run_script.sh R1 BCFTOOLS.nf 10

today=$(date +"%d-%B-%y")
suffix=$1
workflow=$2
threads=$3

for read in ${PWD}/input/*"${suffix}"*
do
	sampleid=$(basename ${read} "${suffix}".fastq.gz)
	echo "nextflow run ${PWD}/${workflow}	--sampleid ${sampleid} --outDir ${PWD}/output/${today} --reads \"${PWD}/input/${sampleid}_{1,2}.fastq.gz\" -resume --threads ${threads}" >> runscript.sh
done

chmod +x runscript.sh

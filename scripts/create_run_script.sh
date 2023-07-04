#!/bin/bash

# bash script to generate a script (run_script_bs.sh) for running in batch mode
# the generated script contains one line for each sample
# the generated script must be run from within the docker container

today=$(date +"%d-%B-%y")

for read in ${PWD}/input/*_1* 
do
	sampleid=$(basename ${read} _1.fastq.gz)
	echo "nextflow run pf.resistance.nf	--sampleid ${sampleid} --outDir ${PWD}/output/${today} --reads \"${PWD}/input/${sampleid}_{1,2}.fastq.gz\" -resume" >> run_script.sh
done

chmod +x run_script_bs.sh

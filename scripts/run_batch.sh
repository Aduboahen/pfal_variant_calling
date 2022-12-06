#!/bin/bash


# USAGE

# bash script running in batch mode without opening the docker container first

# bash run_batch.sh <inputpath> <file_extension> <threads>

# Example:
# bash run_batch.sh /tmp/input/ .fastq.gz 6 

today=$(date +"%d-%B-%y")
thread=${3:-6}

for fname in ${1}/*_1${2}
do
    sample=$(basename ${fname} _1$2 )
    
    docker run -it --rm --mount type=bind,source=${PWD}/output,target=/tmp/output \
    --mount type=bind,source=${1},target=/tmp/input \
   mugicrow/pfal_variant_call:latest \
    nextflow run pf.resistance.nf --sampleid ${sample} \
    --outDir /tmp/output/${today} --threads ${thread} --reads "/tmp/input/${sample}_{1,2}${2}" -resume

done

# pfal_variant_calling
bcftools-based nextflow pipeline to find resistant genes in Plasmodium falciparum WGS

clone this repo

## pull docker image

docker pull mugicrow/pfal_variant_call:latest

copy fastq files into "input" directory

## test 

./script/test_run.sh

## single sample run

./script/start_docker.sh

nextflow run pf.resistance.nf --sampleid test --reads "./input/test_{1,2}.fastq.gz" --outDir ./output -resume

## batch script

./script/run_batch.sh

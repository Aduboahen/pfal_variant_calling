# _P._ _falciparum_ variant calling
bcftools-based nextflow pipeline to find resistant genes in Plasmodium falciparum WGS

## clone this repo

git clone https://github.com/Aduboahen/pfal_variant_call

copy fastq files into "input" directory

## pull docker image

docker pull mugicrow/pfal_variant_call:latest


## test 

./script/test_run.sh

## single sample run

./script/start_docker.sh

nextflow run pf.resistance.nf --sampleid test --reads "./input/test_{1,2}.fastq.gz" --outDir ./output -resume

## batch run

./script/run_batch.sh

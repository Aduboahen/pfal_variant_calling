#!/bin/bash

docker run -i -t --rm \
	--mount type=bind,source=${PWD}/output,target=/tmp \
	mugicrow/pfal_variant_call:latest \
       	nextflow run \
	pf.resistance.nf \
	--sampleid test \
	--outDir ./output \
	--reads "./test_fq/test_{1,2}.fastq.gz" \
	-resume

#!/usr/bin/env python

## python script to generate a script (runscript_py.sh) for running in batch mode
# the generated script contains one line for each sample
# the generated script must be run from within the docker container

import os

from datetime import date

today=date.today().strftime('%Y_%m_%d')

with open("runscript_py.sh", 'w') as outfile:

	for filename in os.listdir(os.getcwd()):
		if filename.endswith("_1.fastq"):
			sample=filename.split('_')
			sampleid=sample[0]
			sample_code=sample[1]
			outfile.write(f'nextflow run fastp.nf --sampleid {sampleid} -resume --outDir ./output --reads "{sampleid}_{{1,2}}.fastq"\n')

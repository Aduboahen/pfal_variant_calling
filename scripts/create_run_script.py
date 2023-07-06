#!/usr/bin/env python

## python script to generate a script (runscript_py.sh) for running in batch mode
# the generated script contains one line for each sample
# the generated script must be run from within the docker container

import os

from datetime import date

today = date.today().strftime('%Y_%m_%d')
pwd = os.getcwd()

with open("runscript.sh", 'w', encoding="utf8") as outfile:
	for filename in os.listdir(os.path.join(pwd, "input")):
		if filename.endswith("_1.fastq.gz"):
			sample = filename.split('_')
			sampleid = sample[0]
			sample_code = sample[1]
			outfile.write(f'nextflow run {os.path.join(pwd, "pf.resistance.nf")} --sampleid {sampleid} --outDir {os.path.join(pwd, "output", today)} --reads "{os.path.join(pwd, "input", sampleid + "_{1,2}.fastq.gz")}" -resume --threads 50 -with-report -with-timeline\n')

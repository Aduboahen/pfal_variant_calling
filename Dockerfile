FROM  mambaorg/micromamba:git-98be205-cuda12.3.2-ubuntu20.04

USER $MAMBA_USER

# install tools
RUN micromamba config append channels bioconda
RUN micromamba config append channels conda-forge
RUN micromamba install -n base -c bioconda -c conda-forge -y python nextflow samtools bcftools bwa picard bedtools mosdepth snpEff fastp ucsc-facount pandas vcfpy pysam multiqc hts-nim-tools

RUN chown $MAMBA_USER:$MAMBA_USER /home/mambauser

# add codebase to docker
COPY --chown=$MAMBA_USER:$MAMBA_USER ./modules/ /home/mambauser/modules

COPY --chown=$MAMBA_USER:$MAMBA_USER ./reference/ /home/mambauser/reference

COPY --chown=$MAMBA_USER:$MAMBA_USER ./ampseq.nf ./ampseq_dsl1.nf ./BCFTOOLS_dsl1.nf ./BCFTOOLS_dsl2.nf ./GATK.nf ./nextflow.config ./scripts/ampseq_batch.sh ./scripts/wgs_batch.sh /home/mambauser/

WORKDIR /home/mambauser

RUN cd /home/mambauser
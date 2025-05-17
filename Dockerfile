FROM  mambaorg/micromamba:git-4427b19-jammy-cuda-11.8.0

USER $MAMBA_USER

# install tools
RUN micromamba config append channels bioconda
RUN micromamba config append channels conda-forge
RUN micromamba install -n base -c bioconda -c conda-forge -y python nextflow samtools bcftools bwa picard bedtools mosdepth snpEff fastp ucsc-facount pandas vcfpy pysam multiqc vim wget

RUN chown $MAMBA_USER:$MAMBA_USER /home/mambauser

# add codebase to docker
COPY --chown=$MAMBA_USER:$MAMBA_USER ./ /home/mambauser

WORKDIR /home/mambauser

FROM  mambaorg/micromamba:git-4427b19-jammy-cuda-11.8.0

USER $MAMBA_USER

# install tools
RUN micromamba install -n base -c bioconda -c conda-forge -y \
 python=3.9 \
 nextflow=22.04.3 \
 samtools \
 bcftools \
 bwa \
 picard \ 
 bedtools \
 mosdepth \
 snpEff \
 fastp \
 ucsc-facount \
 pandas \
 vcfpy \
 pysam \
 multiqc \
 vim \

# add codebase to docker
COPY --chown=$MAMBA_USER:$MAMBA_USER ./ /home/mambauser
WORKDIR /home/mambauser

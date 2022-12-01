FROM ubuntu:22.04

#install curl
RUN apt-get update -y && \
    apt-get upgrade -y && \
    apt-get install -y curl 

# Install miniconda to /miniconda
RUN curl -LO http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -p /miniconda -b && \
    rm Miniconda3-latest-Linux-x86_64.sh

# set the environment
ENV PATH="/miniconda/bin:$PATH" LC_ALL=C

# Install miniconda environment
RUN conda update --all && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda config --add channels anaconda && \
    conda create -n PFAL_RES nextflow samtools bcftools bwa picard bedtools mosdepth snpEff fastp

## activate conda env
RUN echo "source activate PFAL_RES" > ~/.bashrc
ENV PATH=/miniconda/envs/PFAL_RES/bin:${PATH}

# set shell
SHELL ["/bin/bash", "-c"]

# add codebase to docker
ADD ./ /workflow
COPY ./scripts /workflow/scripts
COPY ./pf.resistance.nf /workflow/scripts/pf.resistance.nf
WORKDIR /workflow

FROM continuumio/miniconda3:4.9.2

ENV VERSION "3.1.0"

RUN apt-get --allow-releaseinfo-change update \
    && apt-get install -y git \
    && pip install --upgrade pip \
    && conda update conda \
    && conda install mamba -c conda-forge

RUN apt-get install libarchive13 && mamba install -c bioconda -c conda-forge nanosim=3.1.0 numpy=1.19.5

RUN apt-get install -y procps

WORKDIR /tmp/

ENTRYPOINT ["read_analysis.py"]

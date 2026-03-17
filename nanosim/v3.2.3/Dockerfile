FROM continuumio/miniconda3:4.12.0

ENV VERSION "3.2.3"

RUN apt-get --allow-releaseinfo-change update \
    && apt-get install -y git \
    && pip install --upgrade pip \
    && conda update conda -y \
    && conda install mamba -c conda-forge -y

RUN apt-get install -y libarchive13 && mamba install -c bioconda -c conda-forge nanosim=3.2.3 numpy=1.19.5

RUN apt-get install -y procps

WORKDIR /tmp/

ENTRYPOINT ["read_analysis.py"]

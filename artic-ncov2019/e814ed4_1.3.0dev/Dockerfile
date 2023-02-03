# artic-ncov2019 / fieldbioinformatics
# https://artic.network/ncov-2019/
# Commit: e814ed4f37ec87a194cdc9644cd4b55df27b72de
# https://github.com/artic-network/fieldbioinformatics/commit/00fdbc9d809cb2e33af88fbb14262a1751d49c95
# Commit: 00fdbc9d809cb2e33af88fbb14262a1751d49c95

# base image: miniconda3
FROM continuumio/miniconda3:22.11.1

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

LABEL version="e814ed4f37ec87a194cdc9644cd4b55df27b72de" \
      description="Docker image for artic-ncov2019"

ENV PATH "/opt/conda/envs/artic/bin/:$PATH"

ENV FIELDVERSION "00fdbc9d809cb2e33af88fbb14262a1751d49c95"
ENV NCOVVERSION "e814ed4f37ec87a194cdc9644cd4b55df27b72de"

# Install git for pangolin
RUN apt-get update --fix-missing \
    && apt-get install -y  git build-essential zlib1g zlib1g-dev

#&& python setup.py install \
RUN cd /usr/local/ \
    && conda install -c conda-forge mamba \
    && git clone https://github.com/artic-network/fieldbioinformatics \
    && cd fieldbioinformatics \
    && git checkout $FIELDVERSION \
    && mamba env update  --quiet -f environment.yml \
    && mamba install -n artic -c bioconda -c conda-forge -c defaults artic-network::rampart=1.2.0 snakemake-minimal=5.8.1 \
    && conda init bash && . ~/.bashrc && conda activate artic \
    && pip3 install "setuptools<58.0.0" \
    && python3 -m pip install . \
    && cd ../ \
    && git clone https://github.com/artic-network/artic-ncov2019.git \
    && cd artic-ncov2019 \
    && git checkout $NCOVVERSION

WORKDIR /tmp/

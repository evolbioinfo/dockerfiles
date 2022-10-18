# artic-ncov2019
# https://artic.network/ncov-2019/
# Commit: e814ed4f37ec87a194cdc9644cd4b55df27b72de

# base image: miniconda3
FROM continuumio/miniconda3:4.9.2-alpine

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

LABEL version="e814ed4f37ec87a194cdc9644cd4b55df27b72de" \
      description="Docker image for artic-ncov2019"

ENV PATH "/opt/conda/envs/artic-ncov2019/bin:$PATH"

ENV VERSION "e814ed4f37ec87a194cdc9644cd4b55df27b72de"

# Install git for pangolin
RUN apk update && \
    apk add git bash

RUN cd /usr/local/ \
    && git clone https://github.com/artic-network/artic-ncov2019.git \
    && cd artic-ncov2019 \
    && git checkout $VERSION \
    && conda install -c conda-forge mamba \
    && mamba env remove -n artic-ncov2019 \
    && mamba env create --quiet -f environment.yml \
    && mamba env export --name pangolin > /usr/local/artic-ncov2019/artic.yml 

WORKDIR /tmp/

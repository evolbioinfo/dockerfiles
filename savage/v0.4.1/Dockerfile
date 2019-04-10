# Savage
# https://github.com/haploconduct/haploconduct
# v0.4.1

# base image: miniconda3
FROM continuumio/miniconda3

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>


ENV PATH /opt/conda/envs/env/bin:$PATH

RUN conda update -n base -c defaults conda \
    && conda create -n env python=2.7 \
    && echo "source activate env" > ~/.bashrc \
    && conda config --add channels r \
    && conda config --add channels defaults \
    && conda config --add channels conda-forge \
    && conda config --add channels bioconda \
    && conda install savage=0.4.1

RUN mkdir /pasteur

ENTRYPOINT ["savage"]

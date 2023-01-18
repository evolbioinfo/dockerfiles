# s3cmd
# https://github.com/s3tools/s3cmd
# v2.1.0

# base image: miniconda3
FROM continuumio/miniconda3

# File Author / Maintainer
LABEL maintainer=frederic.lemoine@pasteur.fr


ENV PATH /opt/conda/envs/env/bin:$PATH

RUN conda update -n base -c defaults conda \
    && conda create -n env python=3.5 \
    && echo "source activate env" > ~/.bashrc \
    && pip install s3utils==0.6.1

RUN mkdir /pasteur

ENTRYPOINT ["s3utils"]

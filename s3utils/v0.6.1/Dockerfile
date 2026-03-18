# s3utils
# v0.6.1

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

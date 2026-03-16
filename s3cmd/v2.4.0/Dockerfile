# s3cmd
# https://github.com/s3tools/s3cmd
# v2.4.0

# base image: miniconda3
FROM continuumio/miniconda3

# File Author / Maintainer
LABEL maintainer=frederic.lemoine@pasteur.fr


ENV PATH /opt/conda/envs/env/bin:$PATH

RUN conda update -n base -c defaults conda \
    && conda create -n env python=3.8 \
    && echo "source activate env" > ~/.bashrc \
    && pip install s3cmd==2.4.0 

RUN mkdir /pasteur

ENTRYPOINT ["s3cmd"]

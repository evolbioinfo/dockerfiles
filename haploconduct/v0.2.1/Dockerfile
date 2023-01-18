# HaploConduct
# https://github.com/HaploConduct/HaploConduct/
# v0.2.1

FROM continuumio/miniconda2:4.7.12

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

ENV PATH="/opt/conda/envs/pangolin/bin:$PATH"
ENV VERSION="0.2.1"

#RUN apk update && \
#    apk add git bash
#
RUN cd / \
    && conda config --add channels defaults \
    && conda config --add channels bioconda \
    && conda config --add channels conda-forge \
    && conda install haploconduct=${VERSION} 

WORKDIR /tmp/

ENTRYPOINT ["/opt/conda/bin/haploconduct"]

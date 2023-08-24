# ReSeq
# https://github.com/schmeing/ReSeq
# Commit 053b8d1

# base image: Ubuntu
FROM ubuntu:22.04

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

ENV VERSION=053b8d1

RUN apt-get update --fix-missing \
    && apt-get install -y build-essential zlib1g-dev libbz2-dev python3-dev git cmake libboost-all-dev swig \
    && cd /usr/local/ \
    && git clone https://github.com/schmeing/ReSeq.git \
    && cd ReSeq \
    && git checkout $VERSION \
    && mkdir build \
    && cd build \
    && cmake .. \
    && make \
    && make install \
    && reseq test \
    && cd /usr/local \
    && rm -rf ReSeq \
    && apt-get remove -y git cmake \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["reseq"]

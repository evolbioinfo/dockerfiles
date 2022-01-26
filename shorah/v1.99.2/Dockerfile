# ShoRAH
# https://github.com/cbg-ethz/shorah/
# v1.99.2

# base image: Ubuntu
FROM ubuntu:20.04

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

ENV VERSION=1.99.2
ENV HTSLIBVERSION=1.14
ENV DEBIAN_FRONTEND="noninteractive" TZ="Europe/Paris"
ENV PYTHONPATH="/usr/local/lib/python3.8/site-packages/"

# Intall HTSLIB
RUN apt-get update --fix-missing \
    && apt-get install -y wget gcc make libbz2-dev zlib1g zlib1g-dev liblzma5 liblzma-dev libncurses5 libncurses5-dev bzip2 \
    && cd /usr/local/ \
    && wget https://github.com/samtools/htslib/releases/download/${HTSLIBVERSION}/htslib-${HTSLIBVERSION}.tar.bz2 \ 
    && tar -xjvf htslib-${HTSLIBVERSION}.tar.bz2 \
    && rm -rf htslib-${HTSLIBVERSION}.tar.bz2 \
    && cd htslib-${HTSLIBVERSION} \
    && ./configure \
    && make \
    && make install \
    && cd /usr/local \
    && rm -rf /usr/local/htslib-${HTSLIBVERSION}

# Install SHORAH
RUN apt-get update --fix-missing \
    && apt-get install -y wget gsl-bin libgsl0-dev autoconf automake pkg-config m4 perl zlib1g-dev libboost-all-dev g++ libblas-dev \
    && apt-get install -y python3-biopython python3-numpy \
    && cd /usr/local/ \
    && wget -O shorah-${VERSION}.tar.bz2 https://github.com/cbg-ethz/shorah/releases/download/v${VERSION}/shorah-${VERSION}.tar.bz2 \
    && tar xvjf shorah-${VERSION}.tar.bz2 \
    && rm -f shorah-${VERSION}.tar.bz2 \
    && cd /usr/local/shorah-${VERSION} \
    && ./configure --prefix=/usr/local/ PYTHON=/usr/bin/python3 \
    && make && make install \
    && cd /usr/local \
    && rm -rf shorah-${VERSION} \
    && apt-get remove -y wget libgsl0-dev autoconf automake-1.15 pkg-config \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/shorah"]

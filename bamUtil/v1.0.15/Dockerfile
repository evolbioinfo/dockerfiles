FROM ubuntu:24.04

LABEL maintainer="Julien FUMEY <julien.fumey@pasteur.fr>"

ENV VERSION="v1.0.15"

RUN apt-get update && apt-get install -y \
    clang wget \
    libeigen3-dev \
    libboost-all-dev \
    libomp-dev \
    unzip bzip2 libbz2-dev liblzma-dev \
    cmake \
    g++ \
    git \
    && cd /usr/local/ \
    && git config --global url."https://".insteadOf git:// \
    && git clone --recursive https://github.com/statgen/bamUtil.git \
    && cd bamUtil \
    && git checkout ${VERSION} \
    && make cloneLib \
    && make \
    && cp bin/bam /usr/local/bin \
    && cd /usr/local \
    && wget https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2\
    && tar -xvf samtools-1.21.tar.bz2 \
    && cd samtools-1.21 \
    && ./configure --prefix=/usr/local/ \
    && make \
    && make install \
    && apt-get remove -y wget clang libeigen3-dev libboost-all-dev libomp-dev unzip cmake g++ \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT [ "/bin/bash", "-l", "-c" ]
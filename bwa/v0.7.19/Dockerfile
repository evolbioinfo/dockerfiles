# bwa
# https://github.com/lh3/bwa/
# v0.7.19 (+samtools v1.11)

# base image: Ubuntu
FROM ubuntu:24.04

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

ENV VERSION=0.7.19
ENV SAMVERSION=1.11

RUN apt-get update --fix-missing \
    && apt-get install -y wget gcc make libbz2-dev zlib1g zlib1g-dev liblzma5 liblzma-dev libncurses6 libncurses-dev bzip2 \
    && cd /usr/local/ \
    && wget -O bwa.tar.gz https://github.com/lh3/bwa/archive/refs/tags/v${VERSION}.tar.gz \
    && tar -xzvf bwa.tar.gz \
    && rm -rf bwa.tar.gz \
    && cd bwa-${VERSION} \
    && make \
    && cp bwa *.pl /usr/local/bin/ \
    && cd /usr/local \
    && rm -rf /usr/local/bwa-${VERSION} \
    && wget -O samtools.tar.bz2 https://github.com/samtools/samtools/releases/download/${SAMVERSION}/samtools-${SAMVERSION}.tar.bz2 \
    && tar -xjvf samtools.tar.bz2 \
    && rm -rf samtools.tar.bz2 \
    && cd samtools-${SAMVERSION} \
    && ./configure \
    && make \
    && make install \
    && cd /usr/local \
    && rm -rf /usr/local/samtools-${SAMVERSION} \    
    && apt-get remove -y wget gcc make libbz2-dev zlib1g-dev liblzma-dev libncurses-dev bzip2 \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/bwa"]

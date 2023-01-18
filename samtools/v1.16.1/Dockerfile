# samtools
# http://www.htslib.org/
# v1.16.1

# base image: Ubuntu
FROM ubuntu:22.04

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

ENV VERSION=1.16.1
ENV HTSVERSION=1.16
ENV BCFVERSION=1.16

RUN apt-get update --fix-missing \
    && apt-get install -y wget gcc make libbz2-dev zlib1g zlib1g-dev liblzma5 liblzma-dev libncurses5 libncurses5-dev bzip2 \
    && cd /usr/local/ \
    && wget -O samtools.tar.bz2 https://github.com/samtools/samtools/releases/download/${VERSION}/samtools-${VERSION}.tar.bz2 \
    && tar -xjvf samtools.tar.bz2 \
    && rm -rf samtools.tar.bz2 \
    && cd samtools-${VERSION} \
    && ./configure \
    && make \
    && make install \
    && cd /usr/local \
    && rm -rf /usr/local/samtools-${VERSION} \
    && wget https://github.com/samtools/htslib/releases/download/${HTSVERSION}/htslib-${HTSVERSION}.tar.bz2 \ 
    && tar -xjvf htslib-${HTSVERSION}.tar.bz2 \
    && rm -rf htslib-${HTSVERSION}.tar.bz2 \
    && cd htslib-${HTSVERSION} \
    && ./configure \
    && make \
    && make install \
    && cd /usr/local \
    && rm -rf /usr/local/htslib-${HTSVERSION} \
    && wget https://github.com/samtools/bcftools/releases/download/${BCFVERSION}/bcftools-${BCFVERSION}.tar.bz2 \ 
    && tar -xjvf bcftools-${BCFVERSION}.tar.bz2 \
    && rm -rf bcftools-${BCFVERSION}.tar.bz2 \
    && cd bcftools-${BCFVERSION} \
    && ./configure \
    && make \
    && make install \
    && cd /usr/local \
    && rm -rf /usr/local/bcftools-${BCFVERSION} \
    && apt-get remove -y wget gcc make libbz2-dev zlib1g-dev liblzma-dev libncurses-dev bzip2 \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/samtools"]

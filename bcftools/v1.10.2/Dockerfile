# bcftools
# http://www.htslib.org/
# v1.10.2

# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

ENV VERSION=1.10.2

RUN apt-get update --fix-missing \
    && apt-get install -y wget gcc make libbz2-dev zlib1g zlib1g-dev liblzma5 liblzma-dev libncurses5 libncurses5-dev bzip2 libcurl4-openssl-dev libssl-dev libgetopt-complete-perl  \
    && cd /usr/local/ \
    && wget -O bcftools.tar.bz2 https://github.com/samtools/bcftools/releases/download/${VERSION}/bcftools-${VERSION}.tar.bz2 \
    && tar -xjvf bcftools.tar.bz2 \
    && rm -rf bcftools.tar.bz2 \
    && apt-get install -y libcurl4-openssl-dev libssl-dev libgetopt-complete-perl \
    && wget -O htslib.tar.bz2 https://github.com/samtools/htslib/releases/download/${VERSION}/htslib-${VERSION}.tar.bz2 \
    && tar -xjvf htslib.tar.bz2 \
    && rm htslib.tar.bz2 \
    && cd bcftools-${VERSION} \
    && ./configure \
    && make \
    && make install \
    && rm -rf /usr/local/bcftools-${VERSION} \
    && rm -rf /usr/local/htslib-${VERSION} \
    && apt-get remove -y wget gcc make libbz2-dev zlib1g-dev liblzma-dev libcurl4-openssl-dev libssl-dev libncurses-dev bzip2 \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/bcftools"]

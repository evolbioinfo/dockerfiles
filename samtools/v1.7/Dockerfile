# samtools
# http://www.htslib.org/
# v1.7

# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

RUN apt-get update --fix-missing \
    && apt-get install -y wget gcc make libbz2-dev zlib1g zlib1g-dev liblzma5 liblzma-dev libncurses5 libncurses5-dev bzip2 \
    && cd /usr/local/ \
    && wget -O samtools.tar.bz2 https://github.com/samtools/samtools/releases/download/1.7/samtools-1.7.tar.bz2 \
    && tar -xjvf samtools.tar.bz2 \
    && rm -rf samtools.tar.bz2 \
    && cd samtools-1.7 \
    && ./configure \
    && make \
    && make install \
    && cd /usr/local \
    && rm -rf /usr/local/samtools-1.7 \
    && apt-get remove -y wget gcc make libbz2-dev zlib1g-dev liblzma-dev libncurses-dev bzip2 \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/samtools"]

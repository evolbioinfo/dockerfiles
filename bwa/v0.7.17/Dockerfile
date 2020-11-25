# bwa
# https://github.com/lh3/bwa/
# v0.7.17

# base image: Ubuntu
FROM ubuntu:20.04

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

ENV VERSION=0.7.17


RUN apt-get update --fix-missing \
    && apt-get install -y wget gcc make libbz2-dev zlib1g zlib1g-dev liblzma5 liblzma-dev libncurses5 libncurses5-dev bzip2 \
    && cd /usr/local/ \
    && wget -O bwa.tar.bz2 https://github.com/lh3/bwa/releases/download/v${VERSION}/bwa-${VERSION}.tar.bz2 \
    && tar -xjvf bwa.tar.bz2 \
    && rm -rf bwa.tar.bz2 \
    && cd bwa-${VERSION} \
    && make \
    && cp bwa *.pl /usr/local/bin/ \
    && cd /usr/local \
    && rm -rf /usr/local/bwa-${VERSION} \
    && apt-get remove -y wget gcc make libbz2-dev zlib1g-dev liblzma-dev libncurses-dev bzip2 \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/bwa"]

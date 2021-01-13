# Minimap2
# https://github.com/lh3/minimap2/
# master branch commit 2da649d1d724561d4c2bbe1be9123e2b61bc0029 (16 Nov 2020)

# base image: Ubuntu
FROM ubuntu:20.04

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

ENV VERSION=2da649d1d724561d4c2bbe1be9123e2b61bc0029
ENV SAMTOOLSVERSION=1.11

RUN apt-get update --fix-missing \
    && apt-get install -y wget gcc make libbz2-dev zlib1g zlib1g-dev liblzma5 liblzma-dev libncurses5 libncurses5-dev bzip2 git \
    && cd /usr/local/ \
    && git clone https://github.com/lh3/minimap2.git \
    && cd minimap2 \
    && git checkout ${VERSION} \
    && make \
    && mv minimap2 /usr/local/bin \
    && cd .. \
    && rm -rf minimap2 \
    && cd /usr/local/ \
    && wget -O samtools.tar.bz2 https://github.com/samtools/samtools/releases/download/${SAMTOOLSVERSION}/samtools-${SAMTOOLSVERSION}.tar.bz2 \
    && tar -xjvf samtools.tar.bz2 \
    && rm -rf samtools.tar.bz2 \
    && cd samtools-${SAMTOOLSVERSION} \
    && ./configure \
    && make \
    && make install \
    && cd /usr/local \
    && rm -rf /usr/local/samtools-${SAMTOOLSVERSION} \
    && apt-get remove -y wget gcc make libbz2-dev zlib1g-dev liblzma-dev libncurses-dev bzip2 git \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["minimap2"]

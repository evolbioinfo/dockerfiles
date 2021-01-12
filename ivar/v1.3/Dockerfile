# ivar
# https://github.com/andersen-lab/ivar
# v1.3

# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

ENV SAMTOOLSVERSION=1.11
ENV VERSION=1.3
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib/

RUN apt-get update --fix-missing \
    && apt-get install -y wget gcc make libbz2-dev zlib1g zlib1g-dev liblzma5 liblzma-dev libncurses5 libncurses5-dev bzip2 autoreconf \
    && cd /usr/local/ \
    && wget -O samtools.tar.bz2 https://github.com/samtools/samtools/releases/download/${SAMTOOLSVERSION}/samtools-${SAMTOOLSVERSION}.tar.bz2 \
    && tar -xjvf samtools.tar.bz2 \
    && rm -rf samtools.tar.bz2 \
    && cd samtools-${SAMTOOLSVERSION} \
    && ./configure \
    && make \
    && make install \
    && cd htslib* \
    && ./configure \
    && make \
    && make install \
    && cd /usr/local \
    && rm -rf /usr/local/samtools-${SAMTOOLSVERSION} \
    && wget https://github.com/andersen-lab/ivar/archive/v${VERSION}.tar.gz \
    && tar -xzvf v${VERSION}.tar.gz \
    && cd ivar-${VERSION} \
    && ./autogen.sh && ./configure && make && make install \
    && cd /usr/local/ \
    && rm -rf ivar-${VERSION} \
    && apt-get remove -y wget gcc make libbz2-dev zlib1g-dev liblzma-dev libncurses-dev bzip2 \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/ivar"]

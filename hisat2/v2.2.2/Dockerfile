# HISAT2
# https://github.com/DaehwanKimLab/hisat2
# v2.2.2

# base image: Ubuntu
FROM ubuntu:24.04

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

ENV VERSION=2.2.2
ENV SAMTOOLSVERSION=1.23

RUN apt-get update --fix-missing \
    && apt-get install -y wget g++ gcc make git python3 libbz2-dev zlib1g zlib1g-dev liblzma5 liblzma-dev ncurses-base libncurses-dev bzip2 \
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
    && git clone https://github.com/DaehwanKimLab/hisat2.git \
    && cd hisat2 && git checkout v$VERSION \
    && make \
    && find . -maxdepth 1 -type f -executable -exec sh -c "cp -r {} /usr/local/bin" \; \
    && cp -r scripts/* /usr/local/bin/ \
    && apt-get remove -y wget g++ make git wget gcc g++ make libbz2-dev zlib1g-dev liblzma-dev libncurses-dev bzip2  \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/hisat2"]

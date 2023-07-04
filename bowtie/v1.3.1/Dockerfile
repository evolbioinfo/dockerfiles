# bowtie
# https://github.com/BenLangmead/bowtie/
# v1.3.1

# base image: Ubuntu
FROM ubuntu:22.04

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

ENV VERSION=v1.3.1
ENV SAMTOOLSVERSION=1.17

RUN apt-get update --fix-missing \
    && apt-get install -y git g++ make wget python3 gcc libbz2-dev zlib1g zlib1g-dev liblzma5 liblzma-dev libncurses5 libncurses5-dev bzip2 \
    && cd /usr/local/ \
    && git clone https://github.com/BenLangmead/bowtie.git \
    && cd bowtie \
    && git checkout $VERSION \
    && make && make install \
    && cd /usr/local \
    && rm -rf bowtie \
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
    && apt-get remove -y git gcc make wget libbz2-dev zlib1g-dev liblzma-dev libncurses5-dev bzip2 \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/bowtie"]

# bowtie2
# http://bowtie-bio.sf.net/bowtie2
# v2.5.1

# base image: Ubuntu
FROM ubuntu:22.04

# File Author / Maintainer
LABEL maintainer=frederic.lemoine@pasteur.fr

ENV VERSION=v2.5.1
ENV SAMTOOLSVERSION=1.17

RUN apt-get update --fix-missing \
    && apt-get install -y wget git gcc g++ make libtbb2 libtbb-dev zlib1g zlib1g-dev perl python3 python-is-python3 libbz2-dev zlib1g zlib1g-dev liblzma5 liblzma-dev libncurses5 libncurses5-dev bzip2 \
    && cd /usr/local/ \
    && wget -O bowtie.tar.gz  https://github.com/BenLangmead/bowtie2/archive/${VERSION}.tar.gz \
    && tar -xzvf bowtie.tar.gz \
    && rm -rf bowtie.tar.gz \
    && cd bowtie2-* \
    && make \
    && mv bowtie2 bowtie2-align-s bowtie2-align-l bowtie2-build \
       bowtie2-build-s bowtie2-build-l bowtie2-inspect \
       bowtie2-inspect-s bowtie2-inspect-l /usr/local/bin/ \
    && cd /usr/local \
    && rm -rf /usr/local/bowtie2-* \
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
    && apt-get remove -y wget gcc g++ make libtbb-dev zlib1g-dev libbz2-dev zlib1g-dev liblzma-dev libncurses5-dev bzip2 \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/bowtie2"]

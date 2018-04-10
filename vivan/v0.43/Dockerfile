# ViVan
# http://vivan.tau.ac.il/index.php
# v0.43

# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

# Sources have been downloaded from here:
# http://vivan.tau.ac.il/latest-version.zip
# 2018/04/10
COPY src /usr/local/ViVan

RUN apt-get update --fix-missing \
    && apt-get install -y wget gcc g++ make \
       	       	       	  python2.7 python-pip \
       	       	       	  zlib1g zlib1g-dev \
       	       	       	  libncurses5-dev libncurses5 \
       	       	       	  libbz2-1.0 libbz2-dev \
       	       	       	  liblzma-dev liblzma5 \
    && wget -O bwa.tar.bz2  https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.8.tar.bz2/download \
    && mkdir -p /usr/local/ViVan/bwa/ \
    && tar -xjvf bwa.tar.bz2 \
    && cd bwa-0.7.8/ \
    && make \
    && cp qualfa2fq.pl xa2multi.pl bwa /usr/local/ViVan/bwa/ \
    && cd .. \
    && rm -rf bwa-0.7.8/ bwa.tar.bz2 \
    && pip install biopython scipy numpy  \
    && wget https://github.com/samtools/samtools/releases/download/1.8/samtools-1.8.tar.bz2 \
    && tar -xjvf samtools-1.8.tar.bz2 \
    && cd samtools-1.8 \
    && ./configure --prefix=/usr/local \
    && make && make install \
    && cd .. \
    && rm -rf samtools-1.8 samtools-1.8.tar.bz2 \
    && apt-get remove -y gcc g++ make python-pip zlib1g-dev libncurses5-dev libbz2-dev liblzma-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENV PATH=/usr/local/ViVan/:$PATH

ENTRYPOINT ["/usr/local/ViVan/completeAnalysis.py"]

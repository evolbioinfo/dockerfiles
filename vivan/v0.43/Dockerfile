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
			  libgsl-dev libgsl2 unzip \
    && wget -O bwa.tar.bz2  https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.8.tar.bz2/download \
    && mkdir -p /usr/local/ViVan/bwa/ \
    && tar -xjvf bwa.tar.bz2 \
    && cd bwa-0.7.8/ \
    && make \
    && cp qualfa2fq.pl xa2multi.pl bwa /usr/local/ViVan/bwa/ \
    && cd .. \
    && rm -rf bwa-0.7.8/ bwa.tar.bz2 \
    && pip install biopython scipy numpy  \
    && wget -O samtools-0.1.19.tar.bz2 https://sourceforge.net/projects/samtools/files/samtools/0.1.19/samtools-0.1.19.tar.bz2/download \
    && tar -xjvf samtools-0.1.19.tar.bz2 \
    && cd samtools-0.1.19 \
    && make  \
    && mv samtools /usr/local/bin \
    && cd .. \
    && rm -rf samtools-* \
    && wget -O ea-utils.zip https://github.com/ExpressionAnalysis/ea-utils/archive/bd148d4654fe26d8e218ef8f634fb6a6df9fede1.zip \
    && unzip ea-utils.zip \
    && cd ea-utils-bd148d4654fe26d8e218ef8f634fb6a6df9fede1/clipper \
    && make \
    && mkdir -p /usr/local/ViVan/ea-utils/ \
    && mv multx.sh bwa-to-bowtie determine-phred varcall \
         fastq-multx debug fastq-join fasta-qual-to-fastq \
         fastx-graph contig-stats fastq-to-fasta randomFQ \
         sam-stats qsh merge-fifo varcall-matrix sam-stats.pl \
         gff2gtf grun zhead qsub gtf2bed fastq-clipper \
         check-clipper.sh alc fastq-mcf fastq-stats getgenbankannot \
         /usr/local/ViVan/ea-utils/ \
    && cd ../../ && rm -rf ea-utils* \
    && apt-get remove -y gcc g++ make python-pip zlib1g-dev libncurses5-dev libbz2-dev liblzma-dev libgsl2 unzip \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENV PATH=/usr/local/ViVan/:$PATH

ENTRYPOINT ["/usr/local/ViVan/completeAnalysis.py"]

# bowtie2
# http://bowtie-bio.sf.net/bowtie2
# v2.4.5

# base image: Ubuntu
FROM ubuntu:22.04

# File Author / Maintainer
LABEL maintainer=frederic.lemoine@pasteur.fr

RUN apt-get update --fix-missing \
    && apt-get install -y wget gcc g++ make libtbb2 libtbb-dev zlib1g zlib1g-dev perl python3 python-is-python3 \
    && cd /usr/local/ \
    && wget -O bowtie.tar.gz  https://github.com/BenLangmead/bowtie2/archive/v2.4.5.tar.gz \
    && tar -xzvf bowtie.tar.gz \
    && rm -rf bowtie.tar.gz \
    && cd bowtie2-2.4.5 \
    && make \
    && mv bowtie2 bowtie2-align-s bowtie2-align-l bowtie2-build \
       bowtie2-build-s bowtie2-build-l bowtie2-inspect \
       bowtie2-inspect-s bowtie2-inspect-l /usr/local/bin/ \
    && cd /usr/local \
    && rm -rf /usr/local/bowtie2-2.4.5 \
    && apt-get remove -y wget gcc g++ make libtbb-dev zlib1g-dev \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/bowtie2"]

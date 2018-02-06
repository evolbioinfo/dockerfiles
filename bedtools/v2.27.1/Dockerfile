# bedtools
# http://www.htslib.org/
# 2.27.1

# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

RUN apt-get update --fix-missing \
    && apt-get install -y wget gcc g++ make  zlib1g zlib1g-dev python \
    && cd /usr/local/ \
    && wget -O bedtools.tar.gz https://github.com/arq5x/bedtools2/releases/download/v2.27.1/bedtools-2.27.1.tar.gz \
    && tar -xzvf bedtools.tar.gz \
    && rm -rf bedtools.tar.gz \
    && cd bedtools2 \
    && make \
    && mv bin/* /usr/local/bin \
    && cd /usr/local \
    && rm -rf /usr/local/bedtools2 \
    && apt-get remove -y wget gcc g++ make \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/bedtools"]

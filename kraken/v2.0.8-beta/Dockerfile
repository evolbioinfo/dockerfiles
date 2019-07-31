# Kraken
# https://github.com/DerrickWood/kraken2
# v2.0.8-beta

# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

ENV KRAKEN_VERSION v2.0.8-beta
ENV KRAKEN_DIR=/usr/local/kraken
ENV PATH ${KRAKEN_DIR}:$PATH

RUN apt-get update --fix-missing \
    && apt-get install -y wget make cmake g++ ncbi-blast+=2.2.31-4 \
    && cd /usr/local/ \
    && wget https://github.com/DerrickWood/kraken2/archive/${KRAKEN_VERSION}.tar.gz \
    && tar -xzvf ${KRAKEN_VERSION}.tar.gz \
    && rm -f ${KRAKEN_VERSION}.tar.gz \
    && cd /usr/local/kraken2-*/ \
    && ./install_kraken2.sh ${KRAKEN_DIR} \
    && cd /usr/local/ \
    && rm -rf  kraken2-* \
    && apt-get remove -y make cmake g++ \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ${KRAKEN_DIR}/kraken2

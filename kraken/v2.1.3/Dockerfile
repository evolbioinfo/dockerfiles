# Kraken
# https://github.com/DerrickWood/kraken2
# v2.1.3

# base image: Ubuntu
FROM ubuntu:22.04

# File Author / Maintainer
LABEL maintainer=frederic.lemoine@pasteur.fr

ENV KRAKEN_VERSION v2.1.3
ENV KRAKEN_DIR=/usr/local/kraken
ENV PATH ${KRAKEN_DIR}:$PATH
ENV LC_ALL="en_US.UTF-8"
ENV LANG="en_US.UTF-8"
ENV LANGUAGE="en_US.UTF-8"

RUN apt-get update --fix-missing \
    && apt-get install -y wget make cmake g++ ncbi-blast+=2.12.0+ds-3build1 rsync zlib1g-dev zlib1g \
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

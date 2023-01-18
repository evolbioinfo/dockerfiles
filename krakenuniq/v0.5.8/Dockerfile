# KrakenUniq
# https://github.com/fbreitwieser/krakenuniq
# v0.5.8

# base image: Ubuntu
FROM ubuntu:20.04

# File Author / Maintainer
LABEL maintainer=frederic.lemoine@pasteur.fr

ENV KRAKENUNIQ_VERSION=v0.5.8
ENV JELLYFISH_VERSION=1.1.11
ENV KRAKENUNIQ_PATH=/usr/local/bin/

RUN apt-get update --fix-missing \
    && apt-get install -y wget make g++ gcc git zlib1g zlib1g-dev libgomp1 perl libwww-perl \
    && cd /usr/local/ \
    && wget http://www.cbcb.umd.edu/software/jellyfish/jellyfish-${JELLYFISH_VERSION}.tar.gz \
    && tar xf jellyfish-${JELLYFISH_VERSION}.tar.gz \
    && rm jellyfish-${JELLYFISH_VERSION}.tar.gz \
    && cd jellyfish-${JELLYFISH_VERSION} \
    && ./configure --prefix=/usr/local && make && make install \
    && cd .. && rm -rf jellyfish-${JELLYFISH_VERSION} \
    && git clone https://github.com/fbreitwieser/krakenuniq \
    && cd krakenuniq \
    && git checkout ${KRAKENUNIQ_VERSION} \
    && ./install_krakenuniq.sh ${KRAKENUNIQ_PATH} \
    && cd .. && rm -rf krakenuniq \
    && apt-get remove -y make g++ git \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ${KRAKENUNIQ_PATH}/krakenuniq

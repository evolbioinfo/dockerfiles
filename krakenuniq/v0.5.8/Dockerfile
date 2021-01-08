# KrakenUniq
# https://github.com/fbreitwieser/krakenuniq
# v0.5.8

# base image: Ubuntu
FROM ubuntu:20.04

# File Author / Maintainer
LABEL maintainer=frederic.lemoine@pasteur.fr

ENV KRAKENUNIQ_VERSION=v0.5.8
ENV KRAKENUNIQ_PATH=/usr/local/bin/

RUN apt-get update --fix-missing \
    && apt-get install -y wget make g++ git zlib1g zlib1g-dev \
    && git clone https://github.com/fbreitwieser/krakenuniq \
    && cd krakenuniq \
    && git checkout ${KRAKENUNIQ_VERSION} \
    && ./install_krakenuniq.sh -j ${KRAKENUNIQ_PATH} \
    && cd .. && rm -rf krakenuniq \
    && apt-get remove -y make g++ git \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ${KRAKENUNIQ_PATH}/krakenuniq

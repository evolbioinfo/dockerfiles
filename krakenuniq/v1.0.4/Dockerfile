# KrakenUniq
# https://github.com/fbreitwieser/krakenuniq
# v1.0.4

# base image: Ubuntu
FROM ubuntu:22.04

# File Author / Maintainer
LABEL maintainer=julien.fumey@pasteur.fr

ENV KRAKENUNIQ_VERSION=v1.0.4
ENV KRAKENUNIQ_PATH=/usr/local/bin/

RUN apt-get update --fix-missing \
    && apt-get install -y wget make g++ gcc git zlib1g zlib1g-dev libgomp1 perl libwww-perl libbz2-dev \
    && cd /usr/local/ \
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

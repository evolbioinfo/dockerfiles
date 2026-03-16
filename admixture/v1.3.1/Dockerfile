# admixture
# https://dalexander.github.io/admixture/download.html
# v1.3.1

# base image: ubuntu:20.04
FROM ubuntu:24.04

# File Author / Maintainer
LABEL maintainer frederic.lemoine@pasteur.fr

ENV VERSION=1.3.1

RUN apt-get update --fix-missing \
    && apt-get install -y wget \
    && cd /usr/local \
    && wget https://dalexander.github.io/admixture/binaries/admixture_linux-${VERSION}.tar.gz \
    && tar -xzvf admixture_linux-${VERSION}.tar.gz \
    && mv admixture_linux-${VERSION}/* /usr/local/bin \
    && rm -rf dist admixture_linux-${VERSION}.tar.gz \
    && mkdir /pasteur \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* 

ENTRYPOINT ["admixture"]

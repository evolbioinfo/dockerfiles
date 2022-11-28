# BBMap
# https://sourceforge.net/projects/bbmap/
# v39.01

# base image: Ubuntu
FROM ubuntu:22.04

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

ENV VERSION=39.01

RUN apt-get update --fix-missing \
    && apt-get install -y wget  openjdk-18-jre \
    && cd /usr/local/ \
    && wget -O BBMap_${VERSION}.tar.gz https://sourceforge.net/projects/bbmap/files/BBMap_${VERSION}.tar.gz/download \
    && tar -xzvf BBMap_${VERSION}.tar.gz \
    && rm -rf  BBMap_${VERSION}.tar.gz \
    && mv bbmap/* /usr/local/bin/ \
    && rm -rf bbmap \
    && apt-get remove -y wget gcc \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/bbmap.sh"]

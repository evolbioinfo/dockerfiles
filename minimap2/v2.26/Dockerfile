# Minimap2
# https://github.com/lh3/minimap2/
# v2.26

# base image: Ubuntu
FROM ubuntu:22.04

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

ENV VERSION=2.26

RUN apt-get update --fix-missing \
    && apt-get install -y wget gcc bzip2 make zlib1g zlib1g-dev \
    && cd /usr/local/ \
    && wget https://github.com/lh3/minimap2/releases/download/v${VERSION}/minimap2-${VERSION}.tar.bz2 \
    && tar -xjvf minimap2-${VERSION}.tar.bz2 \
    && cd minimap2-${VERSION}  \
    && make \
    && mv minimap2 /usr/local/bin \
    && cd .. \
    && rm -rf minimap2-${VERSION}.tar.bz2 minimap2-${VERSION} \
    && apt-get remove -y wget gcc bzip2 make zlib1g-dev \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["minimap2"]

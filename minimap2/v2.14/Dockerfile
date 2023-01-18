# Minimap2
# https://github.com/lh3/minimap2/
# v2.14

# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

RUN apt-get update --fix-missing \
    && apt-get install -y wget gcc bzip2 make zlib1g zlib1g-dev \
    && cd /usr/local/ \
    && wget https://github.com/lh3/minimap2/releases/download/v2.14/minimap2-2.14.tar.bz2 \
    && tar -xjvf minimap2-2.14.tar.bz2 \
    && cd minimap2-2.14  \
    && make \
    && mv minimap2 /usr/local/bin \
    && cd .. \
    && rm -rf minimap2-2.14.tar.bz2 minimap2-2.14 \
    && apt-get remove -y wget gcc bzip2 make zlib1g-dev \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["minimap2"]

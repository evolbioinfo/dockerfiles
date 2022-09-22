# MAFFT
# https://mafft.cbrc.jp/alignment/software/
# v7.505

# base image: Ubuntu
FROM ubuntu:22.04

# File Author / Maintainer
LABEL maintainer=frederic.lemoine@pasteur.fr

ENV VERSION=7.505

RUN apt-get update --fix-missing \
    && apt-get install -y wget gcc make \
    && cd /usr/local/ \
    && wget -O mafft-${VERSION}-without-extensions-src.tgz https://mafft.cbrc.jp/alignment/software/mafft-${VERSION}-without-extensions-src.tgz\
    && tar -xzvf mafft-${VERSION}-without-extensions-src.tgz \
    && rm -rf mafft-${VERSION}-without-extensions-src.tgz \
    && cd mafft-${VERSION}-without-extensions/core \
    && make \
    && make install \
    && cd /usr/local \
    && rm -rf /usr/local/mafft-${VERSION}-without-extensions/ \
    && apt-get remove -y wget gcc make \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/mafft"]

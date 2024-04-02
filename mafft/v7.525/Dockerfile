# MAFFT
# https://mafft.cbrc.jp/alignment/software/
# v7.525

# base image: Ubuntu
FROM ubuntu:22.04

# File Author / Maintainer
LABEL maintainer=frederic.lemoine@pasteur.fr

ENV VERSION=7.525

RUN apt-get update --fix-missing \
    && apt-get install -y wget gcc make \
    && cd /usr/local/ \
    && wget -O mafft-${VERSION}-without-extensions-src.tgz https://gitlab.com/sysimm/mafft/-/archive/v${VERSION}/mafft-v${VERSION}.tar.gz \
    && tar -xzvf mafft-${VERSION}-without-extensions-src.tgz \
    && rm -rf mafft-${VERSION}-without-extensions-src.tgz \
    && cd mafft-v${VERSION}/core \
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

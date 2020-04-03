# MAFFT
# https://mafft.cbrc.jp/alignment/software/
# v7.453

# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

RUN apt-get update --fix-missing \
    && apt-get install -y wget gcc make \
    && cd /usr/local/ \
    && wget -O mafft-7.453-without-extensions-src.tgz https://mafft.cbrc.jp/alignment/software/mafft-7.453-without-extensions-src.tgz\
    && tar -xzvf mafft-7.453-without-extensions-src.tgz \
    && rm -rf mafft-7.453-without-extensions-src.tgz \
    && cd mafft-7.453-without-extensions/core \
    && make \
    && make install \
    && cd /usr/local \
    && rm -rf /usr/local/mafft-7.453-without-extensions/ \
    && apt-get remove -y wget gcc make \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/mafft"]

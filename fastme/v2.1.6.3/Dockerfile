# FastME
# https://gite.lirmm.fr/atgc/FastME/
# v2.1.6.3

# base image: Ubuntu
FROM ubuntu:22.04

# File Author / Maintainer
LABEL maintainer=frederic.lemoine@pasteur.fr

ENV VERSION=2.1.6.3

RUN apt-get update --fix-missing \
    && apt-get install -y wget gcc make libgomp1 automake \
    && cd /usr/local/ \
    && wget https://gite.lirmm.fr/atgc/FastME/-/raw/v${VERSION}/tarball/fastme-${VERSION}.tar.gz \
    && tar -xzvf fastme-${VERSION}.tar.gz \
    && cd fastme-${VERSION}/ \
    && ./configure && make && make install \
    && cd .. \
    && rm -rf fastme-${VERSION}.tar.gz \
    && rm -rf fastme* \
    && apt-get remove -y wget g++ make automake \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["fastme"]

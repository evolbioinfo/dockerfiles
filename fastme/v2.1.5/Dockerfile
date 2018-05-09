# FastME
# https://gite.lirmm.fr/atgc/FastME/
# v2.1.5

# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

RUN apt-get update --fix-missing \
    && apt-get install -y wget gcc make libgomp1 automake \
    && cd /usr/local/ \
    && wget http://www.atgc-montpellier.fr/download/sources/fastme/fastme-2.1.5.tar.gz \
    && tar -xzvf fastme-2.1.5.tar.gz \
    && cd fastme-2.1.5 \
    && ./configure && make && make install \
    && cd .. \
    && rm -rf fastme-2.1.5.tar.gz \
    && rm -rf fastme-2.1.5 \
    && apt-get remove -y wget gcc make automake \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["fastme"]

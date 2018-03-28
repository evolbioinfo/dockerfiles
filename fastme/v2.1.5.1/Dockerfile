# FastME
# https://gite.lirmm.fr/atgc/FastME/
# v2.1.5.1

# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

RUN apt-get update --fix-missing \
    && apt-get install -y wget gcc make libgomp1 automake \
    && cd /usr/local/ \
    && wget https://gite.lirmm.fr/atgc/FastME/repository/v2.1.5.1/archive.tar.gz \
    && tar -xzvf archive.tar.gz \
    && cd FastME-v2.1.5.1*/ \
    && ./configure && make && make install \
    && cd .. \
    && rm -rf archive.tar.gz \
    && rm -rf FastME* \
    && apt-get remove -y wget g++ make automake \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["fastme"]

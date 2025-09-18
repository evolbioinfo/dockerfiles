# FastME
# https://gite.lirmm.fr/atgc/FastME/
# v2.1.6.1
# goalign
# https://github.com/evolbioinfo/goalign
# v0.3.1
# gotree
# https://github.com/evolbioinfo/gotree
# v0.2.1

# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

ENV PATH=/usr/local/go/bin:/gopath/bin/:$PATH
ENV GOPATH=/gopath

RUN apt-get update --fix-missing \
    && apt-get install -y wget gcc make libgomp1 automake \
    && cd /usr/local/ \
    && wget https://gite.lirmm.fr/atgc/FastME/-/archive/v2.1.6.1/FastME-v2.1.6.1.tar.gz \
    && tar -xzvf FastME-v2.1.6.1.tar.gz \
    && cd FastME-v2.1.6.1*/ \
    && ./configure && make && make install \
    && cd .. \
    && rm -rf archive.tar.gz \
    && rm -rf FastME*

RUN cd /usr/local/bin \
    && wget https://github.com/evolbioinfo/goalign/releases/download/v0.3.1/goalign_amd64_linux \
    && mv goalign_amd64_linux goalign \
    && chmod +x goalign 

RUN cd /usr/local/bin \
    && wget https://github.com/evolbioinfo/gotree/releases/download/v0.2.10/gotree_amd64_linux \
    && mv gotree_amd64_linux gotree \
    && chmod +x gotree

RUN apt-get remove -y wget g++ make automake \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur


ENTRYPOINT ["fastme"]

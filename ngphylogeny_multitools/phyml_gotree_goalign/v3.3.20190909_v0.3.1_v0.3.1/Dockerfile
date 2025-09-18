# PhyML
# https://github.com/stephaneguindon/phyml/releases
# v3.3.20190909
# goalign
# https://github.com/evolbioinfo/goalign
# v0.3.1
# gotree
# https://github.com/evolbioinfo/gotree
# v0.3.1

# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

ENV PATH=/usr/local/go/bin:/gopath/bin/:$PATH
ENV GOPATH=/gopath

RUN apt-get update --fix-missing \
    && apt-get install -y wget gcc make autoconf automake pkg-config \
    && cd /usr/local/ \
    && wget -O phyml-v3.3.20190909.tar.gz https://github.com/stephaneguindon/phyml/archive/c6916ab672d640225250221be5dded120d910eab.tar.gz \
    && tar -xzvf phyml-v3.3.20190909.tar.gz \
    && rm -rf phyml-v3.3.20190909.tar.gz \
    && cd phyml-c6916ab672d640225250221be5dded120d910eab \
    && sh ./autogen.sh \
    && ./configure \
    && make \
    && make install \
    && cd /usr/local \
    && rm -rf /usr/local/phyml-c6916ab672d640225250221be5dded120d910eab 

RUN cd /usr/local/bin \
    && wget https://github.com/evolbioinfo/goalign/releases/download/v0.3.1/goalign_amd64_linux \
    && mv goalign_amd64_linux goalign \
    && chmod +x goalign 

RUN cd /usr/local/bin \
    && wget https://github.com/evolbioinfo/gotree/releases/download/v0.3.1/gotree_amd64_linux \
    && mv gotree_amd64_linux gotree \
    && chmod +x gotree

RUN apt-get remove -y wget gcc make autoconf automake pkg-config \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/phyml"]

# FastTree
# http://www.microbesonline.org/fasttree/#Install
# v2.1.11
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

ENV VERSION=2.1.11

RUN apt-get update --fix-missing \
    && apt-get install -y wget gcc parallel \
    && cd /usr/local/ \
    && wget http://www.microbesonline.org/fasttree/FastTree-${VERSION}.c \
    && gcc -DNO_SSE -O3 -finline-functions -funroll-loops -Wall -o bin/FastTree FastTree-${VERSION}.c -lm \
    && gcc -DNO_SSE -O3 -finline-functions -funroll-loops -Wall -o bin/FastTree-double -DUSE_DOUBLE FastTree-${VERSION}.c -lm \
    && rm FastTree-${VERSION}.c

RUN cd /usr/local/bin \
    && wget https://github.com/evolbioinfo/goalign/releases/download/v0.3.1/goalign_amd64_linux \
    && mv goalign_amd64_linux goalign \
    && chmod +x goalign 

RUN cd /usr/local/bin \
    && wget https://github.com/evolbioinfo/gotree/releases/download/v0.3.1/gotree_amd64_linux \
    && mv gotree_amd64_linux gotree \
    && chmod +x gotree

RUN apt-get remove -y wget gcc \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/FastTree"]

# FastTree
# http://www.microbesonline.org/fasttree/#Install
# https://github.com/morgannprice/fasttree
# v2.2.0

# base image: Ubuntu
FROM ubuntu:24.04

ENV VERSION="2.2.0"

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

RUN apt-get update --fix-missing \
    && apt-get install -y wget gcc git \
    && cd /usr/local/ \
    && git clone https://github.com/morgannprice/fasttree.git \
    && cd fasttree \
    && git checkout v$VERSION \
    && gcc -DNO_SSE -O3 -finline-functions -funroll-loops -Wall -o /usr/local/bin/FastTree FastTree.c -lm \
    && cd .. && rm -rf fasttree \
    && apt-get remove -y wget gcc git \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/FastTree"]

# LSD2
# https://github.com/tothuhien/lsd2
# v2.3 - https://github.com/tothuhien/lsd2/releases/tag/v2.3

FROM ubuntu:22.04

LABEL maintainer=anna.zhukova@pasteur.fr

ENV VERSION=.2.3

RUN apt-get update --fix-missing \
    && apt-get install -y wget gcc g++ make zip cmake \ 
    && cd /usr/local/ \
    && wget https://github.com/tothuhien/lsd2/archive/v${VERSION}.zip \
    && unzip v${VERSION}.zip \
    && cd lsd2-v${VERSION}/ \
    && cmake . \
    && make \
    && mv lsd2 /usr/local/bin/lsd2 \
    && cd /usr/local \
    && rm -rf v${VERSION}.zip \
    && rm -rf lsd2-v${VERSION}/ \
    && apt-get remove -y wget gcc g++ make zip cmake \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/lsd2"]

# LSD2
# https://github.com/tothuhien/lsd2
# v1.7.1 - https://github.com/tothuhien/lsd2/releases/tag/v1.7.1

FROM ubuntu:16.04

MAINTAINER Anna Zhukova <anna.zhukova@pasteur.fr>

RUN apt-get update --fix-missing \
    && apt-get install -y wget gcc g++ make zip \ 
    && cd /usr/local/ \
    && wget https://github.com/tothuhien/lsd2/archive/v1.7.1.zip \
    && unzip v1.7.1.zip \
    && cd lsd2-1.7.1/src/ \
    && make \
    && mv lsd2 /usr/local/bin/lsd2 \
    && cd /usr/local \
    && rm -rf v1.7.1.zip \
    && rm -rf lsd2-1.7.1/ \
    && apt-get remove -y wget gcc g++ make zip \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/lsd2"]

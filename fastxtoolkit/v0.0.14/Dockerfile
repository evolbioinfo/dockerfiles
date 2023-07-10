# FASTX TOOLKIT
# https://github.com/agordon/fastx_toolkit
# v0.0.14

# base image: Ubuntu
FROM ubuntu:22.04

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

ENV VERSION=0.0.14
COPY patch /usr/local/
COPY patch2 /usr/local/

RUN apt-get update --fix-missing \
    && apt-get install -y wget gcc g++ make libperlio-gzip-perl libgd-graph-perl gnuplot bzip2 patch pkg-config \
    && cd /usr/local/ \
    && wget https://github.com/agordon/libgtextutils/releases/download/0.7/libgtextutils-0.7.tar.gz \
    && tar -xzvf libgtextutils-0.7.tar.gz \
    && cd libgtextutils-0.7 \
    && patch -p1 < ../patch \
    && ./configure && make && make install \
    && cd ../ && rm -rf libgtextutils-0.7* \
    && wget https://github.com/agordon/fastx_toolkit/releases/download/0.0.14/fastx_toolkit-0.0.14.tar.bz2 \
    && tar -xjvf fastx_toolkit-0.0.14.tar.bz2 \
    && cd fastx_toolkit-0.0.14 \
    && patch -p1 < ../patch2 \
    && GTEXTUTILS_LIBS=/usr/local/lib/ && GTEXTUTILS_CFLAGS="-lgtextutils" \
    && ./configure \
    && make  \
    && make install \
    && cd .. && rm -rf fastx_toolkit-0.0.14 \
    && apt-get remove -y wget gcc make autoconf automake g++ \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

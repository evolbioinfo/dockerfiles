# ShoRAH
# https://github.com/cbg-ethz/shorah/
# v1.1.2

# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

RUN apt-get update --fix-missing \
    && apt-get install -y wget gsl-bin libgsl0-dev autoconf automake pkg-config m4 python3 perl zlib1g-dev \
    && apt-get install -y python3-biopython python3-numpy \
    && cd /usr/local/ \
    && wget -O shorah-1.1.2.tar.bz2 https://github.com/cbg-ethz/shorah/releases/download/v1.1.2/shorah-1.1.2.tar.bz2 \
    && tar xvjf shorah-1.1.2.tar.bz2 \
    && rm -f shorah-1.1.2.tar.bz2 \
    && cd /usr/local/shorah-1.1.2 \
    && autoreconf -vif \
    && ./configure \
    && make && make distcheck \
    && make install \
    && apt-get remove -y wget libgsl0-dev autoconf automake-1.15 pkg-config \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* 

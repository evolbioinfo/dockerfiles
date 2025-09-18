# Newick Utilities
# http://cegg.unige.ch/newick_utils
# v1.6
# gotree
# https://github.com/evolbioinfo/gotree
# v0.3.1

# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

RUN apt-get update --fix-missing \
    && apt-get install -y wget gcc make git autoconf flex bison libtool \
    && cd /usr/local/ \
    && git clone https://github.com/tjunier/newick_utils.git \
    && cd newick_utils/ \
    && autoreconf -fi && ./configure && make && make check && make install \
    && cd .. \
    && rm -rf newick* 

RUN cd /usr/local/bin \
    && wget https://github.com/evolbioinfo/gotree/releases/download/v0.3.1/gotree_amd64_linux \
    && mv gotree_amd64_linux gotree \
    && chmod +x gotree

RUN apt-get remove -y wget gcc make git autoconf flex bison libtool \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["nw_display"]

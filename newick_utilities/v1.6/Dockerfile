# Newick Utilities
# http://cegg.unige.ch/newick_utils
# v1.6

# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

RUN apt-get update --fix-missing \
    && apt-get install -y wget gcc make \
    && cd /usr/local/ \
    && wget http://cegg.unige.ch/pub/newick-utils-1.6.tar.gz \
    && tar -xzvf newick-utils-1.6.tar.gz \
    && cd newick-utils-1.6/ \
    && ./configure && make && make install \
    && cd .. \
    && rm -rf newick* \
    && apt-get remove -y wget gcc make \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["nw_display"]

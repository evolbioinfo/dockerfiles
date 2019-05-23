# lastal
# http://last.cbrc.jp/
# v980

# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

RUN apt-get update --fix-missing \
    && apt-get install -y unzip wget make g++ zlib1g zlib1g-dev \
    && wget http://last.cbrc.jp/last-980.zip \
    && unzip last-980.zip \
    && cd last-980/ \
    && make \
    && make install \
    && cd .. \
    && rm -rf last-980.zip last-980 \
    && apt-get remove -y unzip wget make g++ zlib1g-dev \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/lastal"]

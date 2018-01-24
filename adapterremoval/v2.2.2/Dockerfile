# adapter removal
# https://github.com/MikkelSchubert/adapterremoval
# v2.2.2

# base image: Ubuntu
FROM ubuntu:14.04

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

RUN apt-get update --fix-missing \
    && apt-get install -y wget make g++ zlib1g-dev libbz2-dev \
    && cd /usr/local/ \
    && wget -O adapterremoval-2.2.2.tar.gz https://github.com/MikkelSchubert/adapterremoval/archive/v2.2.2.tar.gz \
    && tar xvzf adapterremoval-2.2.2.tar.gz \
    && rm -f adapterremoval-2.2.2.tar.gz \
    && cd /usr/local/adapterremoval-2.2.2 \
    && make && make install \
    && apt-get remove -y wget make g++ zlib1g-dev libbz2-dev \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/AdapterRemoval"]

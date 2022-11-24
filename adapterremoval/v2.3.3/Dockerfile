# adapter removal
# https://github.com/MikkelSchubert/adapterremoval
# v2.3.3

# base image: Ubuntu
FROM ubuntu:22.04

# File Author / Maintainer
LABEL maintainer=frederic.lemoine@pasteur.fr

ENV VERSION=2.3.3

RUN apt-get update --fix-missing \
    && apt-get install -y wget make g++ zlib1g-dev libbz2-dev \
    && cd /usr/local/ \
    && wget -O adapterremoval-${VERSION}.tar.gz https://github.com/MikkelSchubert/adapterremoval/archive/v${VERSION}.tar.gz \
    && tar xvzf adapterremoval-${VERSION}.tar.gz \
    && rm -f adapterremoval-${VERSION}.tar.gz \
    && cd /usr/local/adapterremoval-${VERSION} \
    && make && make install \
    && apt-get remove -y wget make g++ zlib1g-dev libbz2-dev \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/AdapterRemoval"]

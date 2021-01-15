# dsrc
# https://github.com/refresh-bio/DSRC
# v2.0.2

# base image: Ubuntu
FROM ubuntu:20.04

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

ENV VERSION=2.0.2
ENV DEBIAN_FRONTEND="noninteractive" TZ="Europe/Paris"


RUN apt-get update --fix-missing \
    && apt-get install -y wget g++ make libboost-all-dev \
    && cd /usr/local/ \
    && wget https://github.com/refresh-bio/DSRC/archive/v${VERSION}.tar.gz \
    && tar -xzvf v${VERSION}.tar.gz \
    && rm -rf v${VERSION}.tar.gz \
    && cd DSRC-${VERSION} \
    && make bin \
    && mv bin/dsrc /usr/local/bin/ \
    && cd /usr/local \
    && rm -rf /usr/local/DSRC-${VERSION} \
    && apt-get remove -y wget g++ make libboost-all-dev \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/dsrc"]

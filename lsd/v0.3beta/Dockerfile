# LSD
# https://github.com/tothuhien/lsd-0.3beta
# v0.3beta

# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Mathieu Moslonka-Lefebvre <mathieu.moslonka-lefebvre@pasteur.fr>

RUN apt-get update --fix-missing \
    && apt-get install -y subversion gcc g++ make \ 
    && cd /usr/local/ \
    && svn checkout --force https://github.com/tothuhien/lsd-0.3beta/trunk/src \
    && cd src \
    && make \
    && mv lsd /usr/local/bin/lsd \
    && cd /usr/local \
    && rm -rf /usr/local/src \
    && apt-get remove -y subversion gcc g++ make \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/lsd"]

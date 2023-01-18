# LSD
# https://github.com/tothuhien/lsd-0.3beta
# v0.3beta - revision 55183ca9d0ea918314e7c89e3dfd67382b1f2597
# https://github.com/tothuhien/lsd-0.3beta/tree/55183ca9d0ea918314e7c89e3dfd67382b1f2597

# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Mathieu Moslonka-Lefebvre <mathieu.moslonka-lefebvre@pasteur.fr>

RUN apt-get update --fix-missing \
    && apt-get install -y wget gcc g++ make zip \ 
    && cd /usr/local/ \
    && wget https://github.com/tothuhien/lsd-0.3beta/archive/55183ca9d0ea918314e7c89e3dfd67382b1f2597.zip \
    && unzip 55183ca9d0ea918314e7c89e3dfd67382b1f2597.zip \
    && cd lsd-0.3beta-55183ca9d0ea918314e7c89e3dfd67382b1f2597/src \
    && make \
    && mv lsd /usr/local/bin/lsd \
    && cd /usr/local \
    && rm -rf 55183ca9d0ea918314e7c89e3dfd67382b1f2597.zip \
    && rm -rf lsd-0.3beta-55183ca9d0ea918314e7c89e3dfd67382b1f2597/ \
    && apt-get remove -y wget gcc g++ make zip \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/lsd"]

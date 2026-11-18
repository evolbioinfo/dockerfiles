# wget
# https://www.gnu.org/software/wget/
# v1.25.0
# Built from the GNU sources (Ubuntu ships an older wget)

# base image: Ubuntu
FROM ubuntu:24.04

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

ENV DEBIAN_FRONTEND="noninteractive" TZ="Europe/Paris"
ENV VERSION=1.25.0

RUN apt-get update --fix-missing \
    && apt-get install -y curl ca-certificates gcc make pkg-config libssl-dev zlib1g-dev libpcre2-dev libpsl-dev libidn2-dev \
        libssl3t64 zlib1g libpcre2-8-0 libpsl5t64 libidn2-0 \
    && cd /usr/local/src \
    && curl -fsSL -o wget-${VERSION}.tar.gz https://ftp.gnu.org/gnu/wget/wget-${VERSION}.tar.gz \
    && tar -xzf wget-${VERSION}.tar.gz \
    && cd wget-${VERSION} \
    && ./configure --prefix=/usr/local --with-ssl=openssl \
    && make \
    && make install \
    && cd /usr/local/src \
    && rm -rf wget-${VERSION} wget-${VERSION}.tar.gz \
    && apt-get remove -y curl gcc make pkg-config libssl-dev zlib1g-dev libpcre2-dev libpsl-dev libidn2-dev \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/wget"]

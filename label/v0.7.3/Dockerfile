# LABEL
# https://github.com/CDCgov/label
# v0.7.3
# Releases are now published on GitHub (wonder.cdc.gov archives up to v0.6.4)

# base image: Ubuntu
FROM ubuntu:24.04

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

ENV DEBIAN_FRONTEND="noninteractive" TZ="Europe/Paris"
ENV VERSION=0.7.3

RUN apt-get update --fix-missing \
    && apt-get install -y zip unzip perl wget libgomp1 procps \
    && wget https://github.com/CDCgov/label/releases/download/v${VERSION}/label-v${VERSION}-universal.zip \
    && unzip label-v${VERSION}-universal.zip \
    && cd flu-amd/ \
    && mv * /usr/local/bin/ \
    && cd ../ \
    && rm -rf label-v${VERSION}-universal.zip flu-amd \
    && apt-get remove -y wget unzip \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/LABEL"]

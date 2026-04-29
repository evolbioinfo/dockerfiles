# RogueNaRok
# https://github.com/aberer/RogueNaRok
# v1.0.1

# base image: Ubuntu
FROM ubuntu:26.04

# File Author / Maintainer
LABEL maintainer=frederic.lemoine@pasteur.fr

ENV VERSION=v1.0.1
ENV DEBIAN_FRONTEND="noninteractive" TZ="Europe/Paris"
RUN apt-get update --fix-missing \
    && apt-get install -y git make gcc \
    && cd /usr/local/ \
    && git clone --recursive https://github.com/aberer/RogueNaRok.git \
    && cd RogueNaRok \
    && git checkout ${VERSION} \
    && make CFLAGS="-std=gnu17 -Wall -D_GNU_SOURCE -Wl,--allow-multiple-definition" \
    && mv RogueNaRok rnr-lsi rnr-mast rnr-prune rnr-tii /usr/local/bin/ \
    && cd /usr/local/ \
    && rm -rf RogueNaRok \
    && apt-get remove -y git gcc \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/RogueNaRok"]

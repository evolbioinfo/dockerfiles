# seqkit
# https://github.com/shenwei356/seqkit/
# v2.6.1

FROM ubuntu:26.04

ENV VERSION="2.6.1"

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

ENV DEBIAN_FRONTEND="noninteractive" TZ="Europe/Paris"

RUN apt-get update --fix-missing \
    && apt-get install -y wget \
    && wget https://github.com/shenwei356/seqkit/releases/download/v${VERSION}/seqkit_linux_amd64.tar.gz \
    && tar -xzvf seqkit_linux_amd64.tar.gz \
    && rm -f seqkit_linux_amd64.tar.gz \
    && mv seqkit /usr/local/bin/ \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/seqkit"]

# fastqutils
# https://github.com/fredericlemoine/fastqutils
# version v0.1.7

# base image: Ubuntu
FROM ubuntu:24.04

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

ENV PATH=/usr/local/go/bin:/gopath/bin/:$PATH
ENV GOPATH=/gopath
ENV VERSION="0.1.7"
ENV PREFIXVERSION="v"
ENV GOVERSION=1.18
ENV CGO_ENABLED=0

RUN apt-get update --fix-missing \
    && apt-get install -y wget git make xz-utils \
    && cd /usr/local \
    && wget --no-check-certificate https://go.dev/dl/go${GOVERSION}.linux-amd64.tar.gz \
    && tar -C /usr/local -xzf go${GOVERSION}.linux-amd64.tar.gz \
    && rm -f go${GOVERSION}.linux-amd64.tar.gz \
    && mkdir -p /gopath/src/github.com/fredericlemoine/ \
    && cd /gopath/src/github.com/fredericlemoine/ \
    && git clone https://github.com/fredericlemoine/fastqutils.git \
    && cd /gopath/src/github.com/fredericlemoine/fastqutils \
    && git checkout ${PREFIXVERSION}${VERSION} \
    && make  \
    && mv fastqutils /usr/local/bin/ \
    && cd /usr/local \
    && rm -rf /gopath /usr/local/go \
    && apt-get remove -y wget git make \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/fastqutils"]

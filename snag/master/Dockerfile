# snag
# https://github.com/fredericlemoine/snag/
# version master

# base image: Ubuntu
FROM ubuntu:22.04

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

ENV PATH=/usr/local/go/bin:/gopath/bin/:$PATH
ENV GOPATH=/gopath
ENV VERSION=master
ENV PREFIXVERSION=""
ENV GOVERSION=1.21.6
ENV CGO_ENABLED=0

RUN apt-get update --fix-missing \
    && apt-get install -y wget git make \
    && cd /usr/local \
    && wget --no-check-certificate https://storage.googleapis.com/golang/go${GOVERSION}.linux-amd64.tar.gz \
    && tar -C /usr/local -xzf go${GOVERSION}.linux-amd64.tar.gz \
    && rm -f go${GOVERSION}.linux-amd64.tar.gz \
    && mkdir -p /gopath/src/github.com/fredericlemoine/ \
    && cd /gopath/src/github.com/fredericlemoine/ \
    && git clone https://github.com/fredericlemoine/snag.git \
    && cd /gopath/src/github.com/fredericlemoine/snag \
    && git checkout ${PREFIXVERSION}${VERSION} \
    && make \
    && mv snag /usr/local/bin/ \
    && cd /usr/local \
    && rm -rf /gopath /usr/local/go \
    && apt-get remove -y wget git make \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/snag"]

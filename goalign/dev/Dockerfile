# goalign
# https://github.com/evolbioinfo/goalign/
# dev

# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

ENV PATH=/usr/local/go/bin:/gopath/bin/:$PATH
ENV GOPATH=/gopath

RUN apt-get update --fix-missing \
    && apt-get install -y wget git make \
    && cd /usr/local \
    && wget --no-check-certificate https://storage.googleapis.com/golang/go1.10.3.linux-amd64.tar.gz \
    && tar -C /usr/local -xzf go1.10.3.linux-amd64.tar.gz \
    && rm -f go1.10.3.linux-amd64.tar.gz \
    && mkdir -p /gopath/ \
    && go get -u github.com/golang/dep/cmd/dep \
    && mkdir -p /gopath/src/github.com/evolbioinfo/ \
    && cd /gopath/src/github.com/evolbioinfo/ \
    && git clone https://github.com/evolbioinfo/goalign.git \
    && cd /gopath/src/github.com/evolbioinfo/goalign \
    && git checkout dev \
    && make dep && make build \
    && mv goalign /usr/local/bin/ \
    && make test \
    && cd /usr/local \
    && rm -rf /gopath /usr/local/go \
    && apt-get remove -y wget git make \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/goalign"]

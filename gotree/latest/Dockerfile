# gotree
# https://github.com/evolbioinfo/gotree/
# v0.4.1

# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

ENV PATH=/usr/local/go/bin:/gopath/bin/:$PATH
ENV GOPATH=/gopath
ENV VERSION=0.4.1
ENV GOVERSION=1.13.6
ENV CGO_ENABLED=0

RUN apt-get update --fix-missing \
    && apt-get install -y wget git \
    && cd /usr/local \
    && wget --no-check-certificate https://storage.googleapis.com/golang/go${GOVERSION}.linux-amd64.tar.gz \
    && tar -C /usr/local -xzf go${GOVERSION}.linux-amd64.tar.gz \
    && rm -f go${GOVERSION}.linux-amd64.tar.gz \
    && mkdir -p /gopath/ \
    && go get -u github.com/golang/dep/cmd/dep \
    && mkdir -p /gopath/src/github.com/evolbioinfo/ \
    && cd /gopath/src/github.com/evolbioinfo/ \
    && wget -O gotree-${VERSION}.tar.gz https://github.com/evolbioinfo/gotree/archive/v${VERSION}.tar.gz \
    && tar -xzf gotree-${VERSION}.tar.gz \
    && mv gotree-${VERSION} gotree \
    && rm -rf gotree-${VERSION}.tar.gz  \
    && cd /gopath/src/github.com/evolbioinfo/gotree \
    && dep ensure \
    && go build -o gotree -ldflags "-X github.com/evolbioinfo/gotree/cmd.Version=v${VERSION}" github.com/evolbioinfo/gotree \
    && go test github.com/evolbioinfo/gotree/... \
    && bash ./test.sh \
    && mv gotree /usr/local/bin/ \
    && cd /usr/local \
    && rm -rf /gopath /usr/local/go \
    && apt-get remove -y wget git \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/gotree"]

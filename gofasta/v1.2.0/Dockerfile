# gofasta
# https://github.com/virus-evolution/gofasta
# version v1.2.0

# base image: Ubuntu
FROM ubuntu:22.04

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

ENV PATH=/usr/local/go/bin:/gopath/bin/:$PATH
ENV GOPATH=/gopath
ENV VERSION="1.2.0"
ENV PREFIXVERSION="v"
ENV GOVERSION=1.19.6
ENV CGO_ENABLED=0

RUN apt-get update --fix-missing \
    && apt-get install -y wget git xz-utils \
    && cd /usr/local \
    && wget --no-check-certificate https://storage.googleapis.com/golang/go${GOVERSION}.linux-amd64.tar.gz \
    && tar -C /usr/local -xzf go${GOVERSION}.linux-amd64.tar.gz \
    && rm -f go${GOVERSION}.linux-amd64.tar.gz \
    && mkdir -p /gopath/src/github.com/virus-evolution/ \
    && cd /gopath/src/github.com/virus-evolution/ \
    && git clone https://github.com/virus-evolution/gofasta.git \
    && cd /gopath/src/github.com/virus-evolution/gofasta \
    && git checkout ${PREFIXVERSION}${VERSION} \
    && go get . && go build -o gofasta \
    && mv gofasta /usr/local/bin/ \
    && cd /usr/local \
    && rm -rf /gopath /usr/local/go \
    && apt-get remove -y wget git \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/gofasta"]

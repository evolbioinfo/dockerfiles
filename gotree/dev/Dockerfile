# gotree
# https://github.com/evolbioinfo/gotree/
# dev branch

# base image: Ubuntu
FROM ubuntu:20.04

# File Author / Maintainer
LABEL maintainer=frederic.lemoine@pasteur.fr

ENV PATH=/usr/local/go/bin:/gopath/bin/:$PATH
ENV GOPATH=/gopath
ENV VERSION=dev
ENV PREFIXVERSION=""
ENV GOVERSION=1.13.6
ENV CGO_ENABLED=0

RUN apt-get update --fix-missing \
    && apt-get install -y wget git make \
    && cd /usr/local \
    && wget --no-check-certificate https://storage.googleapis.com/golang/go${GOVERSION}.linux-amd64.tar.gz \
    && tar -C /usr/local -xzf go${GOVERSION}.linux-amd64.tar.gz \
    && rm -f go${GOVERSION}.linux-amd64.tar.gz \
    && mkdir -p /gopath/src/github.com/evolbioinfo/ \
    && cd /gopath/src/github.com/evolbioinfo/ \
    && git clone https://github.com/evolbioinfo/gotree.git \
    && cd /gopath/src/github.com/evolbioinfo/gotree \
    && git checkout ${PREFIXVERSION}${VERSION} \
    && make build && make test && make testcommands \
    && mv gotree /usr/local/bin/ \
    && cd /usr/local \
    && rm -rf /gopath /usr/local/go \
    && apt-get remove -y wget git make \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur
    
ENTRYPOINT ["/usr/local/bin/gotree"]

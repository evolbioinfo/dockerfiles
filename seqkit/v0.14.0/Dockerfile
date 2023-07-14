# seqkit
# https://github.com/shenwei356/seqkit/
# v0.14.0

FROM ubuntu:22.04

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

ENV DEBIAN_FRONTEND="noninteractive" TZ="Europe/Paris"

ENV PATH=/usr/local/go/bin:/gopath/bin/:$PATH
ENV GOPATH=/gopath
ENV VERSION="0.14.0"
ENV PREFIXVERSION="v"
ENV GOVERSION=1.20
ENV CGO_ENABLED=0

RUN apt-get update --fix-missing \
    && apt-get install -y wget git make xz-utils \
    && cd /usr/local \
    && wget --no-check-certificate https://storage.googleapis.com/golang/go${GOVERSION}.linux-amd64.tar.gz \
    && tar -C /usr/local -xzf go${GOVERSION}.linux-amd64.tar.gz \
    && rm -f go${GOVERSION}.linux-amd64.tar.gz \
    && mkdir -p /gopath/src/github.com/shenwei356/ \
    && cd /gopath/src/github.com/shenwei356/ \
    && git clone https://github.com/shenwei356/seqkit.git \
    && cd seqkit \
    && git checkout $PREFIXVERSION$VERSION \
    && cd seqkit \
    && go mod init \
    && go mod tidy  \
    && sed -i 's/.*github.com\/shenwei356\/util v0.*//' go.mod \
    && go get github.com/shenwei356/util/byteutil@v0.4.0 \
    && go get github.com/shenwei356/util@v0.4.0  \
    && go build . \
    && mv seqkit /usr/local/bin \
    && cd / && rm -rf gopath /usr/local/go \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/seqkit"]

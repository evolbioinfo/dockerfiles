# gotree
# https://github.com/fredericlemoine/gotree/
# v0.2.7b

# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

ENV PATH=/usr/local/go/bin:/gopath/bin/:$PATH
ENV GOPATH=/gopath

RUN apt-get update --fix-missing \
    && apt-get install -y wget git \
    && cd /usr/local \
    && wget --no-check-certificate https://storage.googleapis.com/golang/go1.9.2.linux-amd64.tar.gz \
    && tar -C /usr/local -xzf go1.9.2.linux-amd64.tar.gz \
    && rm -f go1.9.2.linux-amd64.tar.gz \
    && mkdir -p /gopath/ \
    && go get -u github.com/golang/dep/cmd/dep \
    && mkdir -p /gopath/src/github.com/fredericlemoine/ \
    && cd /gopath/src/github.com/fredericlemoine/ \
    && wget -O gotree-0.2.7b.tar.gz https://github.com/fredericlemoine/gotree/archive/v0.2.7b.tar.gz \
    && tar -xzvf gotree-0.2.7b.tar.gz \
    && mv gotree-0.2.7b gotree \
    && rm -rf gotree-0.2.7b.tar.gz \
    && cd /gopath/src/github.com/fredericlemoine/gotree \
    && dep ensure \
    && go build -o gotree -ldflags "-X github.com/fredericlemoine/gotree/cmd.Version=v0.2.7b" github.com/fredericlemoine/gotree \
    && go test github.com/fredericlemoine/gotree/... \
    && mv gotree /usr/local/bin/ \
    && bash ./test.sh \
    && cd /usr/local \
    && rm -rf /gopath /usr/local/go \
    && apt-get remove -y wget git \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* 

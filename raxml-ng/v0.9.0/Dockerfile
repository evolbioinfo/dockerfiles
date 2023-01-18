# RAxML-ng
# https://github.com/amkozlov/raxml-ng
# v0.9.0

# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

RUN apt-get update --fix-missing \
    && apt-get install -y git flex bison libgmp3-dev cmake make g++ \
    && cd /usr/local/ \
    && git clone --recursive https://github.com/amkozlov/raxml-ng \
    && cd raxml-ng \
    && git checkout 0.9.0 \
    && mkdir build && cd build \
    && cmake .. \
    && make \
    && mv ../bin/raxml-ng /usr/local/bin/ \
    && cd /usr/local/ \
    && rm -rf raxml-ng \
    && apt-get remove -y git flex bison libgmp3-dev make cmake g++ \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/raxml-ng"]

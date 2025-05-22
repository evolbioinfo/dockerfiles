# hyphy
# http://hyphy.org/
# v2.5.73

# base image: Ubuntu
FROM ubuntu:24.04

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

ENV VERSION=2.5.73

RUN apt-get update --fix-missing \
    && apt-get install -y git gcc make cmake libomp-dev g++ zlib1g zlib1g-dev libcurl4t64 libcurl4-openssl-dev libgomp1 \
    && cd /usr/local/ \
    && git clone https://github.com/veg/hyphy.git \
    && cd hyphy \
    && git checkout $VERSION \
    && cmake -DINSTALL_PREFIX=/usr/local/ . \
    && make  \
    && make install \
    && cd /usr/local \
    && rm -rf hyphy \
    && apt-get remove -y git gcc make cmake g++ \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/hyphy"]

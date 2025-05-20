# freebayes
#https://github.com/freebayes/freebayes
# v1.3.9

# base image: Ubuntu
FROM ubuntu:24.04

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

ENV VERSION=1.3.9

RUN apt-get update --fix-missing \
    && apt-get install -y gcc  git bc samtools parallel meson ninja-build libvcflib-dev libvcflib-tools vcftools libwfa2-dev libseqlib-dev cmake \
    && cd /usr/local/ \
    && git clone --recurse-submodules https://github.com/freebayes/freebayes.git \
    && cd freebayes \
    && git checkout v${VERSION} \
    && meson build/ --buildtype release \
    && cd build \
    && ninja \
    && meson test -t 2 \
    && cp freebayes /usr/local/bin \
    && cd /usr/local \
    && rm -rf freebayes \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/freebayes"]

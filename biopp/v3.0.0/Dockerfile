# Bio++
# https://github.com/BioPP/bppsuite
# v3.0.0

# base image: Ubuntu
FROM ubuntu:24.04

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

ENV VERSION=v3.0.0

RUN apt-get update --fix-missing \
    && apt-get install -y git flex bison libgmp3-dev cmake make g++ libeigen3-dev \
    && cd /usr/local/ \
    && for i in bpp-core bpp-seq bpp-phyl bpp-popgen bppsuite; do git clone https://github.com/BioPP/${i} && cd ${i} && git checkout ${VERSION} && cmake -DCMAKE_INSTALL_PREFIX=/usr/local/ ./ &&  make && make install && cd ..; done \
    && apt-get remove -y git flex bison libgmp3-dev make cmake g++ \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/bppml"]

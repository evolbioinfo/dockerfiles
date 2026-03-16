# FastANI
# https://github.com/ParBLiSS/FastANI
# v1.34

# base image: Ubuntu
FROM ubuntu:24.04

# File Author / Maintainer
LABEL maintainer=frederic.lemoine@pasteur.fr

RUN apt-get update --fix-missing \
    && apt-get install -y git autoconf gsl-bin libgsl-dev g++ make zlib1g-dev zlib1g libgomp1 \
    && cd /usr/local/ \
    && git clone https://github.com/ParBLiSS/FastANI \
    && cd FastANI \
    && git checkout v1.34 \
    && sed -i "s/1\.33/1\.34/g" src/map/include/parseCmdArgs.hpp \
    && ./bootstrap.sh \
    && ./configure \
    && make \
    && mv fastANI /usr/local/bin \
    && cd /usr/local \
    && rm -rf FastANI \
    && apt-get remove -y git autoconf libgsl-dev g++ make zlib1g-dev \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/fastANI"]

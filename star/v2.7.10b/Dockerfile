# STAR
# https://github.com/alexdobin/STAR
#  2.7.10b

# base image: Ubuntu
FROM ubuntu:22.04

# File Author / Maintainer
LABEL maintainer=frederic.lemoine@pasteur.fr

ENV VERSION=2.7.10b

RUN apt-get update --fix-missing \
    && apt-get install -y wget make gcc g++ zlib1g zlib1g-dev libgomp1 \
    && cd /usr/local/ \
    && wget https://github.com/alexdobin/STAR/archive/${VERSION}.tar.gz \
    && tar -xzvf ${VERSION}.tar.gz \
    && cd STAR-${VERSION}/source \
    && make STAR \
    && cp STAR /usr/local/bin/ \
    && cd /usr/local/ \
    && rm -rf ${VERSION}.tar.gz STAR-${VERSION} \
    && apt-get remove -y wget zlib1g-dev g++ gcc \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/STAR"]

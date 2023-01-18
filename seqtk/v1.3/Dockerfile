# seqtk
# https://github.com/lh3/seqtk/releases
# v1.3

# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

ENV VERSION=1.3


RUN apt-get update --fix-missing \
    && apt-get install -y wget gcc make libbz2-dev zlib1g zlib1g-dev liblzma5 liblzma-dev libncurses5 libncurses5-dev bzip2 \
    && cd /usr/local/ \
    && wget https://github.com/lh3/seqtk/archive/v${VERSION}.tar.gz \
    && tar -xzvf v${VERSION}.tar.gz \
    && rm v${VERSION}.tar.gz \
    && cd seqtk-${VERSION}/ \
    && make \
    && make install \
    && cd /usr/local \
    && rm -rf /usr/local/seqtk-${VERSION} \
    && apt-get remove -y wget gcc make libbz2-dev zlib1g-dev liblzma-dev libncurses-dev bzip2 \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/seqtk"]

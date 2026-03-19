# Seq-Gen
# https://github.com/rambaut/Seq-Gen/
# http://tree.bio.ed.ac.uk/software/seqgen/
# v1.3.5

# base image: Ubuntu
FROM ubuntu:24.04

ENV VERSION="1.3.5"

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

RUN apt-get update --fix-missing \
    && apt-get install -y wget make gcc \
    && cd /usr/local/ \
    && wget -O Seq-Gen-${VERSION}.tar.gz https://github.com/rambaut/Seq-Gen/archive/refs/tags/v${VERSION}.tar.gz \
    && tar xvzf Seq-Gen-${VERSION}.tar.gz \
    && rm -f Seq-Gen-${VERSION}.tar.gz \
    && cd /usr/local/Seq-Gen-${VERSION}/source \
    && make \
    && mv seq-gen /usr/local/bin/ \
    && cd /usr/local \
    && rm -rf Seq-Gen-${VERSION} \
    && apt-get remove -y wget make gcc \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/seq-gen"]

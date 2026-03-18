# RAxML
# https://github.com/stamatak/standard-RAxML/
# v8.2.13

# base image: Ubuntu
FROM ubuntu:24.04

# File Author / Maintainer
LABEL maintainer=frederic.lemoine@pasteur.fr

ENV VERSION="8.2.13"

RUN apt-get update --fix-missing \
    && apt-get install -y wget gcc make \
    && cd /usr/local/ \
    && wget -O  standard-RAxML-${VERSION}.tar.gz https://github.com/stamatak/standard-RAxML/archive/v${VERSION}.tar.gz \
    && tar -xzvf standard-RAxML-${VERSION}.tar.gz \
    && rm -rf standard-RAxML-${VERSION}.tar.gz \
    && cd standard-RAxML-${VERSION} \
    && make -f Makefile.PTHREADS.gcc \
    && mv raxmlHPC-PTHREADS /usr/local/bin/ \
    && make -f Makefile.AVX2.PTHREADS.gcc \
    && mv raxmlHPC-PTHREADS-AVX2 /usr/local/bin/ \
    && make -f Makefile.SSE3.PTHREADS.gcc \
    && mv raxmlHPC-PTHREADS-SSE3 /usr/local/bin/ \
    && cd /usr/local \
    && rm -rf /usr/local/standard-RAxML-${VERSION} \
    && apt-get remove -y wget gcc make \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/raxmlHPC-PTHREADS"]

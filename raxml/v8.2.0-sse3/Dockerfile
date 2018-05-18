# RAxML
# https://github.com/stamatak/standard-RAxML/
# v8.2.0-PTHREADS-SSE3

# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

RUN apt-get update --fix-missing \
    && apt-get install -y wget gcc make \
    && cd /usr/local/ \
    && wget -O  standard-RAxML-8.2.0.tar.gz https://github.com/stamatak/standard-RAxML/archive/v8.2.0.tar.gz \
    && tar -xzvf standard-RAxML-8.2.0.tar.gz \
    && rm -rf standard-RAxML-8.2.0.tar.gz \
    && cd standard-RAxML-8.2.0 \
    && make -f Makefile.SSE3.PTHREADS.gcc \
    && mv raxmlHPC-PTHREADS-SSE3 /usr/local/bin/raxmlHPC \
    && cd /usr/local \
    && rm -rf /usr/local/standard-RAxML-8.2.0 \
    && apt-get remove -y wget gcc make \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/raxmlHPC"]

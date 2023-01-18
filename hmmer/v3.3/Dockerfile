# hmmer
# http://hmmer.org/
# v3.3

# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

ENV HMMERVERSION="3.3"

ENV PATH="/usr/local/hmmer-3.3/easel/miniapps/:${PATH}"

RUN apt-get update --fix-missing \
    && apt-get install -y wget make gcc \
    && cd /usr/local/ \
    && wget http://eddylab.org/software/hmmer/hmmer-${HMMERVERSION}.tar.gz \
    && tar -xzvf hmmer-${HMMERVERSION}.tar.gz \
    && rm -f hmmer-${HMMERVERSION}.tar.gz \
    && cd hmmer-${HMMERVERSION}/ \
    && ./configure --prefix=/usr/local \
    && make \
    && make install \
    && apt-get remove -y wget make gcc \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/hmmbuild"]

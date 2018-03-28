# Clustal Omega
# http://www.clustal.org/omega/
# v1.2.1

# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

RUN apt-get update --fix-missing \
    && apt-get install -y wget g++ make libargtable2-0 libargtable2-dev file libgomp1 \
    && cd /usr/local/ \
    && wget http://www.clustal.org/omega/clustal-omega-1.2.1.tar.gz \
    && tar -xzvf clustal-omega-1.2.1.tar.gz \
    && rm -f clustal-omega-1.2.1.tar.gz \
    && cd clustal-omega-1.2.1 \
    && ./configure && make && make install \
    && cd .. \
    && rm -rf /usr/local/clustal-omega-1.2.1 \
    && apt-get remove -y wget g++ make libargtable2-dev \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["clustalo"]

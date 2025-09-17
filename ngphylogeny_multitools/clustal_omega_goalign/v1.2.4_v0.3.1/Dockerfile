# Clustal Omega
# http://www.clustal.org/omega/
# v1.2.4
# goalign
# https://github.com/evolbioinfo/goalign
# v0.3.1


# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

ENV PATH=/usr/local/go/bin:/gopath/bin/:$PATH
ENV GOPATH=/gopath

RUN apt-get update --fix-missing \
    && apt-get install -y wget git g++ make libargtable2-0 libargtable2-dev file libgomp1 \
    && cd /usr/local/ \
    && wget http://www.clustal.org/omega/clustal-omega-1.2.4.tar.gz \
    && tar -xzvf clustal-omega-1.2.4.tar.gz \
    && rm -f clustal-omega-1.2.4.tar.gz \
    && cd clustal-omega-1.2.4 \
    && ./configure && make && make install \
    && cd .. \
    && rm -rf /usr/local/clustal-omega-1.2.4

RUN cd /usr/local/bin/ \
    && wget https://github.com/evolbioinfo/goalign/releases/download/v0.3.1/goalign_amd64_linux \
    && mv goalign_amd64_linux goalign \
    && chmod +x goalign \
    && apt-get remove -y wget git g++ make libargtable2-dev \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur


ENTRYPOINT ["clustalo"]

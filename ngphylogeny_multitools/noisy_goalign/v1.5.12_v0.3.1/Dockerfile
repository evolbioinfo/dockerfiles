# Noisy
# http://www.bioinf.uni-leipzig.de/Software/noisy/
# v1.5.12
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
    && apt-get install -y wget git g++ make \
    && cd /usr/local/ \
    && wget http://legacy.bioinf.uni-leipzig.de/Software/noisy/Noisy-1.5.12.tar.gz \
    && tar -xzvf Noisy-1.5.12.tar.gz \
    && rm -rf Noisy-1.5.12.tar.gz \
    && cd Noisy-1.5.12/ \
    && ./configure && make && make install \
    && cd /usr/local \
    && rm -rf Noisy*

RUN cd /usr/local/bin \
    && wget https://github.com/evolbioinfo/goalign/releases/download/v0.3.1/goalign_amd64_linux \
    && mv goalign_amd64_linux goalign \
    && chmod +x goalign \
    && apt-get remove -y wget git g++ make \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur


ENTRYPOINT ["noisy"]

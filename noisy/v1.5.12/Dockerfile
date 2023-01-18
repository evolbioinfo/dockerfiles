# Noisy
# http://www.bioinf.uni-leipzig.de/Software/noisy/
# v1.5.12

# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

RUN apt-get update --fix-missing \
    && apt-get install -y wget g++ make \
    && cd /usr/local/ \
    && wget http://www.bioinf.uni-leipzig.de/Software/noisy/Noisy-1.5.12.tar.gz \
    && tar -xzvf Noisy-1.5.12.tar.gz \
    && rm -rf Noisy-1.5.12.tar.gz \
    && cd Noisy-1.5.12/ \
    && ./configure && make && make install \
    && cd /usr/local \
    && rm -rf Noisy* \
    && apt-get remove -y wget g++ make \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["noisy"]

# Velvet
# https://www.ebi.ac.uk/~zerbino/velvet/
# v1.2.10

# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

ENV PATH="/usr/local/VelvetOptimiser-2.2.4:${PATH}"

RUN apt-get update --fix-missing \
    && apt-get install -y bioperl wget gcc make zlib1g zlib1g-dev \
    && cd /usr/local/ \
    && wget https://www.ebi.ac.uk/~zerbino/velvet/velvet_1.2.10.tgz \
    && tar -xzvf velvet_1.2.10.tgz \
    && cd  velvet_1.2.10/ \
    && make 'MAXKMERLENGTH=249' \
    && mv velvetg velveth /usr/local/bin \
    && mv contrib/VelvetOptimiser-2.2.4 /usr/local/ \
    && cd .. \
    && rm -rf velvet_1.2.10/ \
    && apt-get remove -y wget gcc make zlib1g-dev \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["velveth"]

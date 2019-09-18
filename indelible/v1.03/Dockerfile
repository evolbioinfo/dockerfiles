# INDELIBLE
# http://abacus.gene.ucl.ac.uk/software/indelible/
# v1.03

# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

RUN apt-get update --fix-missing \
    && apt-get install -y git g++ \
    && cd /usr/local/ \
    && git clone https://github.com/evolbioinfo/indelible \
    && cd indelible \
    && git checkout v1.03 \
    && cd src/ \
    && g++ -static -o indelible -m64 -O4 -march=opteron indelible.cpp -lm \
    && mv indelible /usr/local/bin/ \
    && cd /usr/local/ \
    && rm -rf indelible \
    && apt-get remove -y git \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/indelible"]

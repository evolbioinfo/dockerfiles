# PhyML
# http://www.atgc-montpellier.fr/phyml/
# v3.1

# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

RUN apt-get update --fix-missing \
    && apt-get install -y wget unzip \
    && cd /usr/local/ \
    && wget http://www.atgc-montpellier.fr/download/binaries/phyml/PhyML-3.1.zip \
    && unzip PhyML-3.1.zip \
    && rm -rf PhyML-3.1.zip \
    && mv PhyML-3.1/PhyML-3.1_linux64 /usr/local/bin/phyml \
    && cd /usr/local \
    && rm -rf /usr/local/PhyML-3.1 \
    && apt-get remove -y wget unzip \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["phyml"]

# FastCodeML
# https://www2.unil.ch/phylo/bioinformatics/fastcodeml.html
# v1.1.0

# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

RUN apt-get update --fix-missing \
    && apt-get install -y wget \
    && cd /usr/local/ \
    && wget ftp://ftp.vital-it.ch/tools/FastCodeML/FastCodeML-1.1.0.tar.gz \
    && tar -xzvf FastCodeML-1.1.0.tar.gz \
    && rm -rf FastCodeML-1.1.0.tar.gz \
    && mv FastCodeML-1.1.0/fast /usr/local/bin/ \
    && rm -rf FastCodeML-1.1.0 \
    && apt-get remove -y wget \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["fast"]

# BMGE
# https://bioweb.pasteur.fr/packages/pack@BMGE@1.12
# v1.12

# base image: ubuntu:16.04
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>


COPY BMGE.sh /usr/local/bin/BMGE

RUN apt-get update --fix-missing \
    && apt-get install -y wget openjdk-8-jre \
    && cd /usr/local/ \
    && wget ftp://ftp.pasteur.fr/pub/gensoft/projects/BMGE/BMGE-1.12.tar.gz \
    && tar -xzvf BMGE-1.12.tar.gz \
    && rm -rf BMGE-1.12.tar.gz \
    && mv BMGE-1.12/BMGE.jar /usr/local/bin \
    && rm -rf BMGE-1.12 \
    && chmod +x /usr/local/bin/BMGE \
    && apt-get remove -y wget \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/BMGE"]

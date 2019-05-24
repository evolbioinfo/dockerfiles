# IRMA
# https://wonder.cdc.gov/amd/flu/irma/
# v0.6.7

# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

RUN apt-get update --fix-missing \
    && apt-get install -y zip perl wget \
    && wget https://wonder.cdc.gov/amd/flu/irma/flu-amd-201704.zip \
    && unzip flu-amd-201704.zip \
    && cd flu-amd/ \
    && mv * /usr/local/bin/ \
    && cd ../ \
    && rm -rf flu-amd-201704.zip flu-amd \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/IRMA"]

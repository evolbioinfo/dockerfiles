# LABEL
# https://wonder.cdc.gov/amd/flu/label/
# v0.6.4

# base image: r
FROM ubuntu:22.04

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

ENV DEBIAN_FRONTEND="noninteractive" TZ="Europe/Paris"

RUN apt-get update --fix-missing \
    && apt-get install -y zip perl wget libgomp1 \
    && wget https://wonder.cdc.gov/amd/flu/label/flu-amd-LABEL-202209.zip \
    && unzip flu-amd-LABEL-202209.zip \
    && cd flu-amd/ \
    && mv * /usr/local/bin/ \
    && cd ../ \
    && rm -rf  flu-amd-LABEL-202209.zip flu-amd \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/LABEL"]

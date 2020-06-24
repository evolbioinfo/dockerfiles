# NanoPlot
# https://github.com/wdecoster/NanoPlot
# v1.29.1

# base image: Ubuntu
FROM evolbioinfo/python:v3.8.2

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

ENV VERSION=1.29.1

# Switch to your new user in the docker image
USER root

RUN apt-get update \
    && apt-get install -y gcc

# Install treetime
RUN cd /usr/local/ && pip3 install --no-cache-dir NanoPlot==${VERSION} psutil

RUN apt-get remove -y gcc \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* 

# Switch to your new user in the docker image
USER evolbioinfo

# The entrypoint runs pyhton3
ENTRYPOINT ["NanoPlot"]

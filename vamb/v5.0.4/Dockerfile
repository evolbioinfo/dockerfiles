# VAMB
# https://github.com/RasmussenLab/vamb
# v5.0.4

FROM python:3.12-slim

ENV VERSION="5.0.4"

# Python3 needs an UTF-8 locale, http://bugs.python.org/issue19846
ENV LANG C.UTF-8

# Build-time environmental variable so that apt doesn't complain
ARG DEBIAN_FRONTEND=noninteractive

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

# Install build dependencies (required for native extensions: pycoverm, pyrodigal)
RUN apt-get update --fix-missing \
    && apt-get install -y gcc g++ make \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install VAMB
RUN pip3 install --no-cache-dir vamb==${VERSION}

# Clean up
RUN mkdir /pasteur

# The entrypoint runs VAMB with command line arguments
ENTRYPOINT ["vamb"]

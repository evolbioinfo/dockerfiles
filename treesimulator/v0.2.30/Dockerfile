# treesimulator
# https://github.com/evolbioinfo/treesimulator
# v0.2.30
# Based on python:3.12-slim instead of evolbioinfo/python:v3.6.9: treesimulator
# uses namedtuple(defaults=...), which requires python >= 3.7.

FROM python:3.12-slim

ENV VERSION="0.2.30"

# Install bc and treesimulator, create the evolbioinfo user
RUN apt-get update -y \
    && apt-get install -y --no-install-recommends bc \
    && pip3 install --no-cache-dir treesimulator==${VERSION} \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && useradd -m evolbioinfo \
    && mkdir /pasteur

# Switch to your new user in the docker image
USER evolbioinfo

# start working in the "evolbioinfo" home directory
WORKDIR /home/evolbioinfo

# The entrypoint runs command line with command line arguments
ENTRYPOINT ["/bin/bash"]

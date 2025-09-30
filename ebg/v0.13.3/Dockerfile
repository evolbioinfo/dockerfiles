# EBG
# https://github.com/wiegertj/EBG
# v0.13.3
# + RAxML-NG v1.2.2

# base image: Python
FROM python:3.6-slim

# File Author / Maintainer
# From the original Dockerfile of Anna Zhukova PastML
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

# Python3 needs an UTF-8 locale, http://bugs.python.org/issue19846
ENV LANG C.UTF-8
# Build-time environmental variable so that apt doesn't complain
ARG DEBIAN_FRONTEND=noninteractive

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

# Install EBG
RUN cd /usr/local/ && pip3 install --no-cache-dir ebg==0.13.3

# Install RAxML-NG v1.2.2
ENV RAXMLNGVERSION=1.2.2
RUN apt-get update --fix-missing \
    && apt-get install -y git flex bison libgmp3-dev cmake make g++ \
    && cd /usr/local/ \
    && git clone --recursive https://github.com/amkozlov/raxml-ng \
    && cd raxml-ng \
    && git checkout ${RAXMLNGVERSION} \
    && mkdir build && cd build \
    && cmake .. \
    && make \
    && mv ../bin/raxml-ng /usr/local/bin/ \
    && cd /usr/local/ \
    && rm -rf raxml-ng \
    && apt-get remove -y git flex bison libgmp3-dev make cmake g++ \
    && apt-get autoremove -y \
    && apt-get clean

# Clean up
RUN mkdir /pasteur

# create a non-root user named evolbioinfo,
# give them the password "evolbioinfo" put them in the sudo group
RUN useradd -d /home/evolbioinfo -m -s /bin/bash evolbioinfo && echo "evolbioinfo:evolbioinfo" | chpasswd && adduser evolbioinfo sudo

# start working in the "evolbioinfo" home directory
WORKDIR /home/evolbioinfo

# Make the files owned by evolbioinfo
RUN chown -R evolbioinfo:evolbioinfo /home/evolbioinfo

# Switch to your new user in the docker image
USER evolbioinfo

# The entrypoint runs pastml with command line arguments
ENTRYPOINT ["ebg"]

FROM python:3.6-slim

# Python3 needs an UTF-8 locale, http://bugs.python.org/issue19846
ENV LANG C.UTF-8

# Build-time environmental variable so that apt doesn't complain
ARG DEBIAN_FRONTEND=noninteractive

# File Author / Maintainer
MAINTAINER Anna Zhukova <anna.zhukova@pasteur.fr>


RUN apt-get update --fix-missing \
    && apt-get install -y wget gcc g++ make zip

# Install snakemake
RUN cd /usr/local/ && pip3 install --no-cache-dir snakemake==5.6.0

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
ENTRYPOINT ["snakemake"]

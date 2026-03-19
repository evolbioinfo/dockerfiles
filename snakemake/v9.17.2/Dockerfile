FROM python:3.14.3-slim

# Python3 needs an UTF-8 locale, http://bugs.python.org/issue19846
ENV LANG C.UTF-8

ENV VERSION="9.17.2"

# Build-time environmental variable so that apt doesn't complain
ARG DEBIAN_FRONTEND=noninteractive

# File Author / Maintainer
LABEL maintainer=anna.zhukova@pasteur.fr


RUN apt-get update --fix-missing \
    && apt-get install -y wget gcc g++ make zip graphviz libsubid5 \
    && wget https://github.com/apptainer/apptainer/releases/download/v1.5.0-rc.1/apptainer_1.5.0.rc.1-trixie+_amd64.deb \
    && apt-get install -y ./apptainer_1.5.0.rc.1-trixie+_amd64.deb \
    && rm ./apptainer_1.5.0.rc.1-trixie+_amd64.deb \
    && wget https://github.com/apptainer/apptainer/releases/download/v1.5.0-rc.1/apptainer-suid_1.5.0.rc.1-trixie+_amd64.deb \
    && dpkg -i ./apptainer-suid_1.5.0.rc.1-trixie+_amd64.deb \
    && rm apptainer-suid_1.5.0.rc.1-trixie+_amd64.deb

# Install snakemake
RUN cd /usr/local/ && pip3 install --no-cache-dir snakemake==${VERSION}

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

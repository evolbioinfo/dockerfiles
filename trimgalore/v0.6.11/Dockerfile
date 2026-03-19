# TrimGalore
# https://github.com/FelixKrueger/TrimGalore
# v0.6.11

# base image: Ubuntu
FROM ubuntu:24.04

# File Author / Maintainer
LABEL maintainer=frederic.lemoine@pasteur.fr

ENV DEBIAN_FRONTEND="noninteractive" TZ="Europe/Paris"

ENV VERSION="0.6.11"

RUN apt-get update --fix-missing \
    && apt-get install -y wget openjdk-25-jre perl unzip python3-cutadapt cutadapt \
    && cd /usr/local/ \
    && wget -O fastqc.zip https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip \
    && unzip fastqc.zip \
    && rm -rf fastqc.zip \
    && chmod 755 FastQC/fastqc \
    && ln -s /usr/local/FastQC/fastqc /usr/local/bin/fastqc \
    && cd /usr/local \
    && wget https://github.com/FelixKrueger/TrimGalore/archive/${VERSION}.tar.gz \
    && tar -xzvf ${VERSION}.tar.gz \
    && mv TrimGalore-${VERSION}/trim_galore /usr/local/bin \
    && rm -rf ${VERSION}.tar.gz TrimGalore-${VERSION} \
    && apt-get remove -y wget unzip \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/trim_galore"]

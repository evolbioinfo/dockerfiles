# TrimGalore
# https://github.com/FelixKrueger/TrimGalore
# v0.6.6

# base image: Ubuntu
FROM ubuntu:20.04

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>


ENV DEBIAN_FRONTEND="noninteractive" TZ="Europe/Paris"

RUN apt-get update --fix-missing \
    && apt-get install -y wget openjdk-8-jre perl unzip python3-cutadapt cutadapt \
    && cd /usr/local/ \
    && wget -O fastqc.zip http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip \
    && unzip fastqc.zip \
    && rm -rf fastqc.zip \
    && chmod 755 FastQC/fastqc \
    && ln -s /usr/local/FastQC/fastqc /usr/local/bin/fastqc \
    && cd /usr/local \
    && wget https://github.com/FelixKrueger/TrimGalore/archive/0.6.6.tar.gz \
    && tar -xzvf 0.6.6.tar.gz \
    && mv TrimGalore-0.6.6/trim_galore /usr/local/bin \
    && rm -rf 0.6.6.tar.gz TrimGalore-0.6.6 \
    && apt-get remove -y wget unzip \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/trim_galore"]

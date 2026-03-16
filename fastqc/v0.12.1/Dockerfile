# FastQC
# http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
# https://github.com/s-andrews/FastQC/
# v0.12.1

# base image: Ubuntu
FROM ubuntu:24.04

# File Author / Maintainer
LABEL maintainer=frederic.lemoine@pasteur.fr

ENV DEBIAN_FRONTEND="noninteractive" TZ="Europe/Paris"
ENV VERSION=v0.12.1

RUN apt-get update --fix-missing \
    && apt-get install -y wget openjdk-8-jre perl unzip \
    && cd /usr/local/ \
    && wget -O fastqc.zip https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_${VERSION}.zip \
    && unzip fastqc.zip \
    && rm -rf fastqc.zip \
    && chmod 755 FastQC/fastqc \
    && sed -i 's/250/750/g' /usr/local/FastQC/fastqc \
    && ln -s /usr/local/FastQC/fastqc /usr/local/bin/fastqc \
    && apt-get remove -y wget unzip \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/fastqc"]

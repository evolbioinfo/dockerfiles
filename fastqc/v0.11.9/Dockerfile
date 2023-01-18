# FastQC
# http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
# v0.11.9

# base image: Ubuntu
FROM ubuntu:20.04

# File Author / Maintainer
LABEL maintainer=frederic.lemoine@pasteur.fr

ENV DEBIAN_FRONTEND="noninteractive" TZ="Europe/Paris"

RUN apt-get update --fix-missing \
    && apt-get install -y wget openjdk-8-jre perl unzip \
    && cd /usr/local/ \
    && wget -O fastqc.zip http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip \
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

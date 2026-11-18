# FastQC
# http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
# https://github.com/s-andrews/FastQC/
# v0.13.0

# base image: Ubuntu
FROM ubuntu:24.04

# File Author / Maintainer
LABEL maintainer=frederic.lemoine@pasteur.fr

ENV DEBIAN_FRONTEND="noninteractive" TZ="Europe/Paris"
ENV VERSION=v0.13.0

# Since 0.13.0 the fastqc launcher is a python3 script, and the zip is only
# published on GitHub. As in previous images, the default memory per file is
# raised (here from 512MB to 750MB).
RUN apt-get update --fix-missing \
    && apt-get install -y wget openjdk-8-jre python3 unzip \
    && cd /usr/local/ \
    && wget -O fastqc.zip https://github.com/s-andrews/FastQC/releases/download/${VERSION}/fastqc_${VERSION}.zip \
    && unzip fastqc.zip \
    && rm -rf fastqc.zip \
    && chmod 755 FastQC/fastqc \
    && sed -i 's/Defaults to 512MB/Defaults to 750MB/; s/default=512)/default=750)/' /usr/local/FastQC/fastqc \
    && grep -q 'default=750)' /usr/local/FastQC/fastqc \
    && ln -s /usr/local/FastQC/fastqc /usr/local/bin/fastqc \
    && apt-get remove -y wget unzip \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/fastqc"]

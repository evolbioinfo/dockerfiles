# TrimGalore
# https://github.com/FelixKrueger/TrimGalore
# v2.5.0
# Note: Trim Galore 2.x is a single binary (Rust rewrite) with Cutadapt and
# FastQC built in: they are no longer installed in the image.
# See https://www.trimgalore.com/reference/migration/

# base image: Ubuntu
FROM ubuntu:24.04

# File Author / Maintainer
LABEL maintainer=frederic.lemoine@pasteur.fr

ENV DEBIAN_FRONTEND="noninteractive" TZ="Europe/Paris"
ENV VERSION="2.5.0"

RUN apt-get update --fix-missing \
    && apt-get install -y wget \
    && cd /usr/local \
    && wget https://github.com/FelixKrueger/TrimGalore/releases/download/v${VERSION}/trim_galore-linux-x86_64.tar.gz \
    && tar -xzvf trim_galore-linux-x86_64.tar.gz \
    && mv trim_galore-linux-x86_64/trim_galore /usr/local/bin/ \
    && chmod +x /usr/local/bin/trim_galore \
    && rm -rf trim_galore-linux-x86_64 trim_galore-linux-x86_64.tar.gz \
    && apt-get remove -y wget \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/trim_galore"]

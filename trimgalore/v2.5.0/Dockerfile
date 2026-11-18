# TrimGalore
# https://github.com/FelixKrueger/TrimGalore
# v2.5.0
# Note: Trim Galore 2.x is a single binary (Rust rewrite) with Cutadapt and
# FastQC built in: they are no longer installed in the image.
# See https://www.trimgalore.com/reference/migration/
# The binary needs glibc >= 2.43, hence Ubuntu 26.04. It is downloaded in an
# Ubuntu 24.04 stage: tar in Ubuntu 26.04 fails under Rosetta amd64 emulation
# when building on Apple Silicon.

FROM ubuntu:24.04 AS download

ENV VERSION="2.5.0"

RUN apt-get update --fix-missing \
    && apt-get install -y wget \
    && cd /tmp \
    && wget https://github.com/FelixKrueger/TrimGalore/releases/download/v${VERSION}/trim_galore-linux-x86_64.tar.gz \
    && tar -xzvf trim_galore-linux-x86_64.tar.gz \
    && chmod +x trim_galore-linux-x86_64/trim_galore

# base image: Ubuntu
FROM ubuntu:26.04

# File Author / Maintainer
LABEL maintainer=frederic.lemoine@pasteur.fr

ENV DEBIAN_FRONTEND="noninteractive" TZ="Europe/Paris"
ENV VERSION="2.5.0"

COPY --from=download /tmp/trim_galore-linux-x86_64/trim_galore /usr/local/bin/trim_galore

RUN mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/trim_galore"]

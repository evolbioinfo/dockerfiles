# RNA-Seqc
# https://software.broadinstitute.org/cancer/cga/rna-seqc
# Archived at https://github.com/francois-a/rnaseqc/
# v1.1.9

FROM ubuntu:22.04

ENV VERSION="1.1.9"

COPY rnaseqc.sh /usr/local/bin
RUN apt-get update --fix-missing \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y wget openjdk-11-jre unzip \
    && cd /usr/local/bin/ \
    && wget https://github.com/francois-a/rnaseqc/releases/download/v1.1.9/RNA-SeQC_1.1.9.zip \
    && unzip RNA-SeQC_1.1.9.zip \
    && rm RNA-SeQC_1.1.9.zip \
    && chmod +x /usr/local/bin/rnaseqc.sh \
    && apt-get remove -y wget unzip \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/rnaseqc.sh"]
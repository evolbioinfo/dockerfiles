# RNA-Seqc
# https://software.broadinstitute.org/cancer/cga/rnaseqc_download
# v1.1.8

FROM ubuntu:22.04

ENV VERSION="1.1.8"

COPY rnaseqc.sh /usr/local/bin
RUN apt-get update --fix-missing \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y wget openjdk-11-jre unzip \
    && cd /usr/local/bin/ \
    && wget http://www.broadinstitute.org/cancer/cga/tools/rnaseqc/RNA-SeQC_v${VERSION}.jar \
    && chmod +x /usr/local/bin/rnaseqc.sh \
    && apt-get remove -y wget unzip \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/rnaseqc.sh"]
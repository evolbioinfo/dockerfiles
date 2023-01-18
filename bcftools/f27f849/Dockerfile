# bcftools
# http://www.htslib.org/
# commit f27f849

# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

ENV VERSION=f27f849
ENV HTSVERSION=1.10.2

RUN apt-get update --fix-missing \
    && apt-get install -y wget gcc make git libbz2-dev zlib1g zlib1g-dev liblzma5 liblzma-dev libncurses5 libncurses5-dev bzip2 libcurl4-openssl-dev libssl-dev libgetopt-complete-perl  \
    && cd /usr/local/ \
    && git clone https://github.com/samtools/htslib.git \
    && cd htslib && git checkout ${HTSVERSION} \
    && cd .. \
    && git clone https://github.com/samtools/bcftools.git \
    && cd bcftools \
    && git checkout ${VERSION} \
    && make \
    && make install \
    && rm -rf /usr/local/bcftools \
    && rm -rf /usr/local/htslib \
    && apt-get remove -y wget gcc make libbz2-dev zlib1g-dev liblzma-dev libncurses-dev bzip2 \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/bcftools"]

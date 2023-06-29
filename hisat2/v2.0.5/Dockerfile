# HISAT2
# https://github.com/DaehwanKimLab/hisat2
# v2.0.5

# base image: Ubuntu
FROM ubuntu:22.04

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

ENV VERSION=2.0.5

COPY patch /usr/local/

RUN apt-get update --fix-missing \
    && apt-get install -y wget g++ make git patch python2 \
    && ln -s /usr/bin/python2 /usr/bin/python \
    && cd /usr/local/ \
    && git clone https://github.com/DaehwanKimLab/hisat2.git \
    && cd hisat2 && git checkout v$VERSION \
    && patch < ../patch \
    && make \
    && find . -maxdepth 1 -type f -executable -exec sh -c "cp -r {} /usr/local/bin" \; \
    && cp -r scripts/* /usr/local/bin/ \
    && apt-get remove -y wget g++ make git patch \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/hisat2"]

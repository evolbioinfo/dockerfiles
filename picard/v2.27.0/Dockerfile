# Picard Tools
# https://github.com/broadinstitute/picard
# v2.27.0

# base image: ubuntu:20.04
FROM ubuntu:20.04

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

COPY picard /usr/local/bin/picard

RUN apt-get update --fix-missing \
    && apt-get install -y wget openjdk-8-jre \
    && cd /usr/local/bin/ \
    && wget https://github.com/broadinstitute/picard/releases/download/2.27.0/picard.jar \
    && chmod +x /usr/local/bin/picard \
    && apt-get remove -y wget \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/picard"]

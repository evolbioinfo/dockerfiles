# Picard Tools
# https://github.com/broadinstitute/picard
# v3.4.0

# base image: ubuntu:24.04
FROM ubuntu:24.04

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"
ENV VERSION=3.4.0
COPY picard /usr/local/bin/picard

RUN apt-get update --fix-missing \
    && apt-get install -y wget openjdk-25-jre \
    && cd /usr/local/bin/ \
    && wget https://github.com/broadinstitute/picard/releases/download/${VERSION}/picard.jar \
    && chmod +x /usr/local/bin/picard \
    && apt-get remove -y wget \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/picard"]

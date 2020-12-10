# DamageProfiler
# https://github.com/Integrative-Transcriptomics/DamageProfiler
# v1.1

# base image: ubuntu:20.04
FROM ubuntu:20.04

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

COPY DamageProfiler.sh /usr/local/bin/DamageProfiler

ENV JAVA_OPTS="-Xmx2g"

RUN apt-get update --fix-missing \
    && apt-get install -y wget openjdk-14-jre less\
    && cd /usr/local/ \
    && wget https://github.com/Integrative-Transcriptomics/DamageProfiler/releases/download/v1.1/DamageProfiler-1.1-java11.jar \
    && chmod +x /usr/local/bin/DamageProfiler \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/DamageProfiler"]

# jq
# https://jqlang.org/
# v1.8.2

# base image: ubuntu:16.04
FROM ubuntu:24.04

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

RUN apt-get update --fix-missing \
    && apt-get install -y wget curl \
    && wget -O /usr/local/bin/jq https://github.com/jqlang/jq/releases/download/jq-1.8.2/jq-linux64 \
    && chmod +x /usr/local/bin/jq \
    && mkdir /pasteur \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* 

ENTRYPOINT ["/usr/local/bin/jq"]

# Alfred
# https://github.com/tobiasrausch/alfred
# v0.6.2

# base image: Ubuntu
FROM ubuntu:24.04

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

ENV ALFRED_VERSION=v0.6.2

# Static binary published by upstream
RUN apt-get update --fix-missing \
    && apt-get install -y wget \
    && wget -O /usr/local/bin/alfred https://github.com/tobiasrausch/alfred/releases/download/${ALFRED_VERSION}/alfred-${ALFRED_VERSION}-linux-amd64 \
    && chmod +x /usr/local/bin/alfred \
    && apt-get remove -y wget \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/alfred"]

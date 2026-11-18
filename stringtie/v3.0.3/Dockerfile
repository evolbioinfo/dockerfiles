# stringtie
# https://github.com/gpertea/stringtie
# v3.0.3

# base image: Ubuntu
FROM ubuntu:24.04

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

ENV VERSION=3.0.3

# Linux binary published by upstream (prepDE.py3 needs python3)
RUN apt-get update --fix-missing \
    && apt-get install -y wget python3 \
    && cd /usr/local/ \
    && wget https://github.com/gpertea/stringtie/releases/download/v${VERSION}/stringtie-${VERSION}.Linux_x86_64.tar.gz \
    && tar -xzvf stringtie-${VERSION}.Linux_x86_64.tar.gz \
    && mv stringtie-${VERSION}.Linux_x86_64/stringtie stringtie-${VERSION}.Linux_x86_64/prepDE.py3 /usr/local/bin/ \
    && chmod +x /usr/local/bin/stringtie /usr/local/bin/prepDE.py3 \
    && rm -rf stringtie-${VERSION}.Linux_x86_64 stringtie-${VERSION}.Linux_x86_64.tar.gz \
    && apt-get remove -y wget \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/stringtie"]

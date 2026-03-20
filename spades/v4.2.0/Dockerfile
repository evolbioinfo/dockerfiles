# SPAdes
# http://cab.spbu.ru/software/spades/
# https://github.com/ablab/spades
# v4.2.0

# base image: Ubuntu
FROM ubuntu:24.04

# File Author / Maintainer
LABEL maintainer=frederic.lemoine@pasteur.fr

ENV VERSION="4.2.0"

RUN apt-get update --fix-missing \
    && apt-get install -y wget libbz2-1.0 zlib1g bzip2 python3 python-is-python3 libgomp1 \
    && cd /usr/local/ \
    && wget https://github.com/ablab/spades/releases/download/v${VERSION}/SPAdes-${VERSION}-Linux.tar.gz \
    && tar -xzf SPAdes-${VERSION}-Linux.tar.gz \
    && rm -rf SPAdes-${VERSION}-Linux.tar.gz \
    && mv SPAdes-${VERSION}-Linux/bin/* /usr/local/bin/ \
    && mv SPAdes-${VERSION}-Linux/share/* /usr/local/share/ \
    && rm -rf SPAdes-${VERSION}-Linux \
    && apt-get remove -y  g++ cmake libbz2-dev zlib1g-dev \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/spades.py"]

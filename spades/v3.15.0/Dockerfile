# SPAdes
# http://cab.spbu.ru/software/spades/
# v3.15.0

# base image: Ubuntu
FROM ubuntu:20.04

# File Author / Maintainer
LABEL maintainer=frederic.lemoine@pasteur.fr

ENV VERSION="3.15.0"

RUN apt-get update --fix-missing \
    && apt-get install -y wget g++ cmake libbz2-1.0 libbz2-dev zlib1g zlib1g-dev bzip2 python2.7 python  libgomp1 \
    && cd /usr/local/ \
    && wget http://cab.spbu.ru/files/release${VERSION}/SPAdes-${VERSION}.tar.gz \
    && tar -xzf SPAdes-${VERSION}.tar.gz \
    && cd SPAdes-${VERSION} \
    && PREFIX=/usr/local ./spades_compile.sh \
    && cd .. \
    && rm -rf SPAdes-${VERSION} \
    && apt-get remove -y  g++ cmake libbz2-dev zlib1g-dev \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/spades.py"]

# SPAdes
# http://cab.spbu.ru/software/spades/
# v3.10.0

# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

RUN apt-get update --fix-missing \
    && apt-get install -y wget g++ cmake libbz2-1.0 libbz2-dev zlib1g zlib1g-dev bzip2 python2.7 python  libgomp1 \
    && cd /usr/local/ \
    && wget http://cab.spbu.ru/files/release3.10.0/SPAdes-3.10.0.tar.gz \
    && tar -xzf SPAdes-3.10.0.tar.gz \
    && cd SPAdes-3.10.0 \
    && PREFIX=/usr/local ./spades_compile.sh \
    && cd .. \
    && rm -rf SPAdes-3.10.0 \
    && apt-get remove -y  g++ cmake libbz2-dev zlib1g-dev \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/spades.py"]

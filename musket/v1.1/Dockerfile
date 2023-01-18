# musket
# http://musket.sourceforge.net/homepage.htm
# v1.1

# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

RUN apt-get update --fix-missing \
    && apt-get install -y wget g++ make zlib1g zlib1g-dev libgomp1 \
    && cd /usr/local/ \
    && wget -O musket-1.1.tar.gz https://sourceforge.net/projects/musket/files/musket-1.1.tar.gz/download \
    && tar -xzvf musket-1.1.tar.gz \
    && rm -rf musket-1.1.tar.gz \
    && cd musket-1.1 \
    && make \
    && cp musket /usr/local/bin \
    && cd /usr/local \
    && rm -rf /usr/local/musket-1.1 \
    && apt-get remove -y wget g++ make zlib1g-dev \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/musket"]

# Khmer 
# http://khmer.readthedocs.io
# v2.0

# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

RUN apt-get update --fix-missing \
    && apt-get install -y python2.7-dev python-virtualenv python-pip gcc g++ \
    && pip install khmer==2.0 \
    && apt-get remove -y gcc g++ \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/normalize-by-median.py"]

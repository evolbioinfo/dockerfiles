# Seq-Gen
# https://github.com/rambaut/Seq-Gen/
# http://tree.bio.ed.ac.uk/software/seqgen/
# v1.3.4

# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

RUN apt-get update --fix-missing \
    && apt-get install -y wget make gcc \
    && cd /usr/local/ \
    && wget -O Seq-Gen-1.3.4.tar.gz https://github.com/rambaut/Seq-Gen/archive/1.3.4.tar.gz \
    && tar xvzf Seq-Gen-1.3.4.tar.gz \
    && rm -f Seq-Gen-1.3.4.tar.gz \
    && cd /usr/local/Seq-Gen-1.3.4/source \
    && make \
    && mv seq-gen /usr/local/bin/ \
    && cd /usr/local \
    && rm -rf Seq-Gen-1.3.4 \
    && apt-get remove -y wget make gcc \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* 

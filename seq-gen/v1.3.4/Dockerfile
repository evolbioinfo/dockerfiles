# Seq-Gen
# https://github.com/rambaut/Seq-Gen/
# http://tree.bio.ed.ac.uk/software/seqgen/
# v1.3.4

# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>
RUN apt-get update --fix-missing
RUN apt-get install -y wget make gcc

WORKDIR /usr/local/
RUN wget -O Seq-Gen-1.3.4.tar.gz https://github.com/rambaut/Seq-Gen/archive/1.3.4.tar.gz \
    && tar xvzf Seq-Gen-1.3.4.tar.gz \
    && rm -f Seq-Gen-1.3.4.tar.gz
WORKDIR /usr/local/Seq-Gen-1.3.4/source
RUN make
RUN mv seq-gen /usr/local/bin/
WORKDIR /usr/local
RUN rm -rf Seq-Gen-1.3.4

RUN apt-get remove -y wget make gcc

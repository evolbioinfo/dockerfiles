# FastTree
# http://www.microbesonline.org/fasttree/#Install
# v2.1.10

# base image: Ubuntu
FROM ubuntu:14.04

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

RUN apt-get update --fix-missing
RUN apt-get install -y wget gcc

WORKDIR /usr/local/
RUN wget http://www.microbesonline.org/fasttree/FastTree-2.1.10.c
RUN gcc -DNO_SSE -O3 -finline-functions -funroll-loops -Wall -o bin/FastTree FastTree-2.1.10.c -lm
RUN rm FastTree-2.1.10.c

RUN apt-get remove -y wget gcc \
    && apt-get clean

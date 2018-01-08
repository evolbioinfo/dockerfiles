# adapter removal
# https://github.com/MikkelSchubert/adapterremoval
# v2.2.2

# base image: Ubuntu
FROM ubuntu:14.04

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

RUN apt-get update --fix-missing
RUN apt-get install -y wget make g++ zlib1g-dev libbz2-dev

WORKDIR /usr/local/
RUN wget -O adapterremoval-2.2.2.tar.gz https://github.com/MikkelSchubert/adapterremoval/archive/v2.2.2.tar.gz \
    && tar xvzf adapterremoval-2.2.2.tar.gz \
    && rm -f adapterremoval-2.2.2.tar.gz
WORKDIR /usr/local/adapterremoval-2.2.2
RUN make && make install

RUN apt-get remove -y wget make g++ zlib1g-dev libbz2-dev \
    && apt-get clean

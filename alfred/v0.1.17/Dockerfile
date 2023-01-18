# Alfred
# https://github.com/tobiasrausch/alfred
# v0.1.17

# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

ENV BOOST_ROOT /usr
ENV ALFRED_VERSION v0.1.17

# install required packages
RUN apt-get update \
    && apt-get install -y build-essential cmake g++ gfortran git \
        hdf5-tools libboost-date-time-dev libboost-program-options-dev libboost-system-dev \
	libboost-filesystem-dev libboost-iostreams-dev libbz2-dev libhdf5-dev libncurses-dev \
	liblzma-dev zlib1g-dev \
    && cd /usr/local \
    && git clone --recursive https://github.com/tobiasrausch/alfred.git \
    && cd alfred \
    && git checkout ${ALFRED_VERSION} \
    && make STATIC=1 all \
    && cp src/alfred /usr/local/bin \
    && cd .. \
    && rm -rf alfred \
    && apt-get remove -y build-essential cmake g++ gfortran git libboost-date-time-dev \
       libboost-program-options-dev libboost-system-dev libboost-filesystem-dev \
       libboost-iostreams-dev libbz2-dev libhdf5-dev libncurses-dev liblzma-dev zlib1g-dev \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/alfred"]

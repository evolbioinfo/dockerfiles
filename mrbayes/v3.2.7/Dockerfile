# mrbayes
# https://github.com/NBISweden/MrBayes
# v3.2.7

# base image: Ubuntu
FROM ubuntu:20.04

# File Author / Maintainer
LABEL maintainer=frederic.lemoine@pasteur.fr
ENV DEBIAN_FRONTEND="noninteractive" TZ="Europe/Paris"
RUN apt-get update --fix-missing \
    && apt-get install -y wget gcc g++ make automake autoconf libtool ocl-icd-opencl-dev openjdk-11-jdk \
    && cd /usr/local/ \
    && wget https://github.com/beagle-dev/beagle-lib/archive/refs/tags/v3.1.2.tar.gz \
    && tar -xzvf v3.1.2.tar.gz \
    && cd beagle-lib-3.1.2 \
    && ./autogen.sh && ./configure && make && make install \
    && cd ../ \
    && rm -rf beagle* \
    && wget https://github.com/NBISweden/MrBayes/releases/download/v3.2.7/mrbayes-3.2.7.tar.gz \
    && tar -xzvf mrbayes-3.2.7.tar.gz \
    && cd mrbayes-3.2.7/ \
    && ./configure && make && make install \
    && cd /usr/local \
    && rm -rf mrbayes* \
    && apt-get remove -y wget gcc g++ make automake autoconf \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["mb"]

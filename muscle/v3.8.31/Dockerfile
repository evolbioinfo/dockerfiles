# MUSCLE
# http://www.drive5.com/muscle/downloads.htm
# v3.8.31

# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

RUN apt-get update --fix-missing \
    && apt-get install -y wget g++ make \
    && cd /usr/local/ \
    && wget http://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_src.tar.gz \
    && tar -xzvf muscle3.8.31_src.tar.gz \
    && cd muscle3.8.31/src \
    && make \
    && mv muscle /usr/local/bin \
    && cd ../../ \
    && rm -rf muscle* \
    && apt-get remove -y wget g++ make \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["muscle"]

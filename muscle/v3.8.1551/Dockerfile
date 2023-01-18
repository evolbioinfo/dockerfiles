# MUSCLE
# http://www.drive5.com/muscle/downloads.htm
# v3.8.1551

# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

RUN apt-get update --fix-missing \
    && apt-get install -y wget g++ make \
    && mkdir /usr/local/muscle/ \
    && cd /usr/local/muscle/ \
    && wget http://www.drive5.com/muscle/muscle_src_3.8.1551.tar.gz \
    && tar -xzvf muscle_src_3.8.1551.tar.gz \
    && make \
    && mv muscle /usr/local/bin \
    && cd .. \
    && rm -rf muscle \
    && apt-get remove -y wget g++ make \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["muscle"]

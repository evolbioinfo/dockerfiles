# Trimal
# https://github.com/scapella/trimal/
# v1.4.1

# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

RUN apt-get update --fix-missing \
    && apt-get install -y wget g++ make  \
    && cd /usr/local/ \
    && wget https://github.com/scapella/trimal/archive/v1.4.1.tar.gz \
    && tar -xzvf v1.4.1.tar.gz \
    && cd trimal-1.4.1/source \
    && make \
    && mv readal statal trimal /usr/local/bin \
    && cd ../.. \
    && rm -rf  v1.4.1.tar.gz trimal-1.4.1  \
    && apt-get remove -y wget g++ make \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["trimal"]

# STAR
# https://github.com/alexdobin/STAR
# 2.7.6a

# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

RUN apt-get update --fix-missing \
    && apt-get install -y wget \
    && cd /usr/local/ \
    && wget https://github.com/alexdobin/STAR/archive/2.7.6a.tar.gz \
    && tar -xzvf 2.7.6a.tar.gz \
    && cp STAR-2.7.6a/bin/Linux_x86_64_static/* /usr/local/bin \
    && rm -rf 2.7.6a.tar.gz STAR-2.7.6a \
    && apt-get remove -y wget \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/STAR"]

# canu
# http://canu.readthedocs.io/
# 2.3

# base image: Ubuntu
FROM ubuntu:24.04

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

ENV VERSION v2.3
ENV PATH=/usr/local/bin:$PATH

RUN apt-get update --fix-missing \
    && apt-get install -y libfindbin-libs-perl wget xz-utils \
    && wget https://github.com/marbl/canu/releases/download/v2.3/canu-2.3.Linux-amd64.tar.xz \
    && tar --directory /usr/local --strip-components=1 -xJvf canu-2.3.Linux-amd64.tar.xz canu-2.3/ \ 
    && rm canu-2.3.Linux-amd64.tar.xz \
    && apt-get remove -y wget \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["canu"]

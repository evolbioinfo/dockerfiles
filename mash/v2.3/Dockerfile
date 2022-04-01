# Mash
# http://mash.readthedocs.io/en/latest/
# v2.3

# base image: Ubuntu
FROM ubuntu:20.04

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

RUN apt-get update --fix-missing \
    && apt-get install -y wget capnproto autoconf libgsl23 libgsl-dev g++ zlib1g zlib1g-dev libcapnp-0.7.0 libcapnp-dev make \
    && ln -s /usr/lib/x86_64-linux-gnu/libcapnp.a /usr/lib/libcapnp.a \
    && ln -s /usr/lib/x86_64-linux-gnu/libkj.a /usr/lib/libkj.a \
    && cd /usr/local/ \
    && wget https://github.com/marbl/Mash/archive/v2.3.tar.gz \
    && tar -xzf  v2.3.tar.gz \
    && cd Mash-2.3 \
    && ./bootstrap.sh \
    && ./configure --prefix=/usr/local --with-capnp=/usr \
    && make \
    && make install \
    && cd .. \
    && rm -rf Mash-2.3 \
    && apt-get remove -y gcc libgsl-dev zlib1g-dev autoconf libcapnp-dev \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/mash"]

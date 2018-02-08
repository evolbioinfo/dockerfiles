# PhyML
# https://github.com/stephaneguindon/phyml/releases
# v3.3.20180129

# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

RUN apt-get update --fix-missing \
    && apt-get install -y wget gcc make autoconf automake pkg-config \
    && cd /usr/local/ \
    && wget -O phyml-3.3.20180129.tar.gz https://github.com/stephaneguindon/phyml/archive/3.3.20180129.tar.gz \
    && tar -xzvf phyml-3.3.20180129.tar.gz \
    && rm -rf phyml-3.3.20180129.tar.gz \
    && cd phyml-3.3.20180129 \
    && sh ./autogen.sh \
    && ./configure \
    && make \
    && make install \
    && cd /usr/local \
    && rm -rf /usr/local/phyml-3.3.20180129 \
    && apt-get remove -y wget gcc make autoconf automake pkg-config \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/phyml"]

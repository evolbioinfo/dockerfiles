# PhyML
# https://github.com/stephaneguindon/phyml/releases

# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

ENV PHYML_VERSION=v3.3.20190321
ENV PHYML_COMMIT=fd7e21ec129611e0750c4654f21edcfecf25e7f7

RUN apt-get update --fix-missing \
    && apt-get install -y wget gcc make autoconf automake pkg-config \
    && cd /usr/local/ \
    && wget -O phyml-${PHYML_VERSION}.tar.gz https://github.com/stephaneguindon/phyml/archive/${PHYML_COMMIT}.tar.gz \
    && tar -xzvf phyml-${PHYML_VERSION}.tar.gz \
    && rm -rf phyml-${PHYML_VERSION}.tar.gz \
    && cd phyml-${PHYML_COMMIT} \
    && sh ./autogen.sh \
    && ./configure \
    && make \
    && make install \
    && cd /usr/local \
    && rm -rf /usr/local/phyml-${PHYML_COMMIT} \
    && apt-get remove -y wget gcc make autoconf automake pkg-config \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/phyml"]

# PhyML
# https://github.com/stephaneguindon/phyml/releases
# v3.3.20260528

# base image: Ubuntu
FROM ubuntu:26.04

# File Author / Maintainer
LABEL maintainer=frederic.lemoine@pasteur.fr

ENV PHYMLVERSION=v3.3.20260528
ENV PHYMLCOMMIT=42bb4b52f479938a3427dfe4357b861f0e306686
ENV DEBIAN_FRONTEND="noninteractive" TZ="Europe/Paris"

RUN apt-get update --fix-missing \
    && apt-get install -y wget gcc make autoconf automake pkg-config \
    && cd /usr/local/ \
    && wget -O phyml-${PHYMLVERSION}.tar.gz https://github.com/stephaneguindon/phyml/archive/${PHYMLCOMMIT}.tar.gz \
    && tar -xzvf phyml-${PHYMLVERSION}.tar.gz \
    && rm -rf phyml-${PHYMLVERSION}.tar.gz \
    && cd phyml-${PHYMLCOMMIT} \
    && sh ./autogen.sh \
    && ./configure \
    && make \
    && make install \
    && cd /usr/local \
    && rm -rf /usr/local/phyml-${PHYMLCOMMIT} \
    && apt-get remove -y wget gcc make autoconf automake pkg-config \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/phyml"]

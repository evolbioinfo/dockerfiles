# lastal
# https://gitlab.com/mcfrith/last
# v1654

# base image: Ubuntu
FROM ubuntu:24.04

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

ENV VERSION=1654

# last.cbrc.jp no longer hosts the archives: sources come from GitLab.
# python3 is needed by the helper scripts installed with LAST.
RUN apt-get update --fix-missing \
    && apt-get install -y wget make g++ zlib1g zlib1g-dev python3 \
    && cd /usr/local/ \
    && wget https://gitlab.com/mcfrith/last/-/archive/${VERSION}/last-${VERSION}.tar.gz \
    && tar -xzvf last-${VERSION}.tar.gz \
    && cd last-${VERSION}/ \
    && make \
    && make install \
    && cd .. \
    && rm -rf last-${VERSION}.tar.gz last-${VERSION} \
    && apt-get remove -y wget make g++ zlib1g-dev \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/lastal"]

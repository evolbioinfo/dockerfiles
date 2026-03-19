# Trimal
# https://github.com/scapella/trimal/
# v1.5.1

# base image: Ubuntu
FROM ubuntu:24.04

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

ENV VERSION="1.5.1"

RUN apt-get update --fix-missing \
    && apt-get install -y wget g++ make  \
    && cd /usr/local/ \
    && wget https://github.com/scapella/trimal/archive/v${VERSION}.tar.gz \
    && tar -xzvf v${VERSION}.tar.gz \
    && cd trimal-${VERSION}/source \
    && make \
    && mv readal statal trimal /usr/local/bin \
    && cd ../.. \
    && rm -rf  v${VERSION}.tar.gz trimal-${VERSION}  \
    && apt-get remove -y wget g++ make \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["trimal"]

# epa-ng
# https://github.com/Pbdas/epa-ng
# v0.3.8

# base image: Ubuntu
FROM ubuntu:20.04

# File Author / Maintainer
LABEL maintainer=adrien.le-meur@pasteur.fr

ENV DEBIAN_FRONTEND="noninteractive" TZ="Europe/Paris"

RUN apt-get update --fix-missing \
    && apt-get install -y wget g++ make autotools-dev libtool flex bison cmake automake autoconf zlib1g-dev libgomp1 \
    && export CXX="g++" \
    && cd /usr/local/ \
    && wget https://github.com/Pbdas/epa-ng/archive/v0.3.8.tar.gz \
    && tar -xzvf v0.3.8.tar.gz \
    && rm v0.3.8.tar.gz \
    && cd epa-ng-0.3.8 \
    && make \
    && cp bin/* /usr/local/bin/ \
    && rm -rf /usr/local/epa-ng-0.3.8 \
    && apt-get remove -y wget g++ make autotools-dev libtool flex bison cmake automake autoconf \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/epa-ng"]

# canu
# http://canu.readthedocs.io/
# 2.0

# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

ENV VERSION v2.0
ENV PATH=/usr/local/bin:$PATH

RUN apt-get update --fix-missing \
    && apt-get install -y git make g++ libfindbin-libs-perl \
    && git clone https://github.com/marbl/canu.git \
    && cd canu \
    && git checkout ${VERSION} \
    && cd src \
    && make \
    && cd ../ \
    && cp -r Linux-amd64/* /usr/local/ \
    && cd .. \
    && rm -rf canu \
    && apt-get remove -y git make g++ \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["canu"]

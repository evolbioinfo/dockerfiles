# fineStructure
# https://people.maths.bris.ac.uk/~madjl/finestructure/finestructure.html
# v4.1.1

# base image: Ubuntu
FROM ubuntu:20.04

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

ENV VERSION=4.1.1
ENV CHROMOVERSION=0.0.4
ENV PATH="/usr/local/fs_${VERSION}/:${PATH}"

RUN apt-get update --fix-missing \
    && apt-get install -y wget make unzip libswitch-perl libgomp1 libgsl23 libgsl-dev g++ \
    && cd /usr/local/ \
    && wget https://people.maths.bris.ac.uk/~madjl/finestructure/fs_${VERSION}.zip \
    && unzip fs_${VERSION}.zip \
    && rm fs_${VERSION}.zip \
    && mv fs_${VERSION}/fs_linux* fs_${VERSION}/fs \
    && wget https://people.maths.bris.ac.uk/~madjl/finestructure-old/chromocombine-${CHROMOVERSION}.tar.gz \
    && tar -xzvf chromocombine-${CHROMOVERSION}.tar.gz \
    && cd chromocombine-${CHROMOVERSION} \
    && ./configure \
    && make && make install \
    && cd ../ \
    && rm -rf chromocombine-${CHROMOVERSION} \
    && rm chromocombine-${CHROMOVERSION}.tar.gz \
    && apt-get remove -y wget unzip make libgsl-dev g++ \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/fs_4.1.1/fs"]

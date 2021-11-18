# fineStructure
# https://people.maths.bris.ac.uk/~madjl/finestructure/finestructure.html
# v4.1.1

# base image: Ubuntu
FROM ubuntu:20.04

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

ENV VERSION=4.1.1
ENV CHROMOVERSION=0.0.4
ENV CHROMOPAINTERVERSION=9eb996ed5e7a9c5aa2c339edbf812deb48ac81af
ENV PATH="/usr/local/fs_${VERSION}/:${PATH}"
ENV DEBIAN_FRONTEND="noninteractive" TZ="Europe/Paris"
ENV CRANREPO="'http://cran.univ-paris1.fr/'"
ENV RLIBPATH="'/usr/local/lib/R/site-library/'"

RUN apt-get update --fix-missing \
    && apt-get install -y git r-base wget make unzip libswitch-perl libgomp1 libgsl23 libgsl-dev g++ zlib1g-dev zlib1g \
    && R -e "install.packages('optparse',lib=$RLIBPATH,repo=$CRANREPO)"  \
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
    && git clone https://github.com/sahwa/ChromoPainterV2.git \
    && cd ChromoPainterV2 \
    && git checkout 9eb996ed5e7a9c5aa2c339edbf812deb48ac81af \
    && gcc -Wall -o ChromoPainterv2 ChromoPainterv2.c -lm -lz \
    && mv ChromoPainterv2 impute2chromopainter2.pl /usr/local/bin \
    && cd .. && rm -rf ChromoPainterV2 \
    && git clone https://github.com/sahwa/sourcefindV2.git \
    && cd sourcefindV2 \
    && git checkout a823fe71427c95175637bb6db2e4ff3019d0a079 \
    && chmod +x sourcefindv2.R && cp sourcefindv2.R /usr/local/bin/ \
    && cd ../ && rm -rf sourcefindV2 \
    && apt-get remove -y libgsl-dev g++ zlib1g-dev git \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/fs_4.1.1/fs"]

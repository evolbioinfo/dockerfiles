# Schmutzi
# https://github.com/grenaud/schmutzi/
# v1.5.6

# base image: Ubuntu
FROM r-base:4.1.1

# File Author / Maintainer
LABEL maintainer=frederic.lemoine@pasteur.fr

ENV CRANREPO="'https://pbil.univ-lyon1.fr/CRAN/'"
ENV RFORGEREPO="'http://R-Forge.R-project.org'"
ENV RLIBPATH="'/usr/local/lib/R/site-library/'"
ENV SCHMUTZIVERSION=v1.5.6
ENV DEBIAN_FRONTEND="noninteractive" TZ="Europe/Paris"
RUN apt-get update --fix-missing \
    && apt-get install -y perl git cmake g++-8 zlib1g zlib1g-dev \
    && cd /usr/local/ \
    && R -e "install.packages('MASS', lib=$RLIBPATH, repos=$CRANREPO)" \
    && R -e "install.packages('fitdistrplus', lib=$RLIBPATH, repos=$CRANREPO)" \
    && git clone https://github.com/grenaud/schmutzi.git \
    && cd schmutzi \
    && git clone https://github.com/pezmaster31/bamtools.git lib/bamtools && cd lib/bamtools && git checkout v2.5.1 && cd ../../ \
    && git clone https://github.com/grenaud/libgab.git lib/libgab && cd lib/libgab && git checkout 1.0.4 && cd ../../ \
    && git checkout ${SCHMUTZIVERSION} \
    && make && make install \
    && cd .. \
    && rm -rf schmutzi \
    && apt-get remove -y git cmake g++-8 zlib1g-dev \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/bin/schmutzi.pl"]

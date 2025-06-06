# BayesTraits
# http://www.evolution.reading.ac.uk/BayesTraitsV4.0.0/BayesTraitsV4.0.0.html
# v4.0.0

# base image: Ubuntu
FROM ubuntu:24.04

# File Author / Maintainer
LABEL maintainer=frederic.lemoine@pasteur.fr

RUN apt-get update --fix-missing \
    && apt-get install -y wget \
    && cd /usr/local/ \
    && wget http://www.evolution.reading.ac.uk/BayesTraitsV5.0.0/Files/BayesTraitsV5.0.0-Linux.tar.gz \
    && tar -xzvf BayesTraitsV5.0.0-Linux.tar.gz \
    && rm -rf BayesTraitsV5.0.0-Linux.tar.gz \
    && cp BayesTraitsV5.0.0-Linux/BayesTraitsV5 /usr/local/bin/ \
    && rm -rf BayesTraitsV5.0.0-Linux \
    && apt-get remove -y wget \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["BayesTraitsV5"]

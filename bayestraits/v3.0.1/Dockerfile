# BayesTraits
# http://www.evolution.reading.ac.uk/BayesTraitsV3.0.1/BayesTraitsV3.0.1.html
# v3.0.1

# base image: Ubuntu
FROM ubuntu:22.04

# File Author / Maintainer
LABEL maintainer=frederic.lemoine@pasteur.fr

RUN apt-get update --fix-missing \
    && apt-get install -y wget \
    && cd /usr/local/ \
    && wget http://www.evolution.reading.ac.uk/BayesTraitsV3.0.1/Files/BayesTraitsV3.0.1-Linux.tar.gz \
    && tar -xzvf BayesTraitsV3.0.1-Linux.tar.gz \
    && rm -rf BayesTraitsV3.0.1-Linux.tar.gz \
    && cp BayesTraitsV3.0.1-Linux/BayesTraitsV3 /usr/local/bin/ \
    && rm -rf BayesTraitsV3.0.1-Linux \
    && apt-get remove -y wget \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["BayesTraitsV3"]

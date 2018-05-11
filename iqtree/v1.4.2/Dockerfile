# IQ-TREE
# http://www.iqtree.org/
# v1.4.2

# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

RUN apt-get update --fix-missing \
    && apt-get install -y wget libgomp1 gcc cmake make eigen3 libeigen3-dev \
    && cd /usr/local/ \
    && wget https://github.com/Cibiv/IQ-TREE/archive/v1.4.2.tar.gz \
    && tar -xzvf v1.4.2.tar.gz \
    && cd IQ-TREE-1.4.2 \
    && mkdir build && cd build \
    && cmake -DIQTREE_FLAGS=omp .. \
    && make \
    && mv iqtree-omp /usr/local/bin/iqtree \
    && ../../ \
    && rm -rf v1.4.2.tar.gz \
    && rm -rf IQ-TREE-1.4.2 \
    && apt-get remove -y wget cmake make gcc libeigen3-dev \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["iqtree"]

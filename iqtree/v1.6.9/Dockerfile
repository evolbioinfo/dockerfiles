# IQ-TREE
# http://www.iqtree.org/
# v1.6.9

# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

RUN apt-get update --fix-missing \
    && apt-get install -y wget libgomp1 gcc cmake make libeigen3-dev \
    && cd /usr/local/ \
    && wget https://github.com/Cibiv/IQ-TREE/archive/v1.6.9.tar.gz \
    && tar -xzvf v1.6.9.tar.gz \
    && cd IQ-TREE-1.6.9 \
    && mkdir build && cd build \
    && cmake -DIQTREE_FLAGS=omp .. \
    && make \
    && mv iqtree /usr/local/bin/iqtree \
    && cd ../../ \
    && rm -rf v1.6.9.tar.gz \
    && rm -rf IQ-TREE-1.6.9 \
    && apt-get remove -y wget cmake make gcc \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["iqtree"]

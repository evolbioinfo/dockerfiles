# GUBBINS
# https://github.com/nickjcroucher/gubbins
# v3.4.3

FROM ubuntu:26.04

LABEL maintainer="Julien FUMEY <julien.fumey@pasteur.fr>"

ENV RAXML_VERSION="8.2.13"
ENV IQTREE_VERSION="3.1.0"
ENV RAXMLNG_VERSION="2.0.0"
ENV FASTTREE_VERSION="2.2.0"
ENV RAPIDNJ=""

RUN apt-get update --fix-missing && apt-get install -y \
    python3 python3-numba python3-biopython python3-dendropy python3-scipy python3-numpy python3-multiprocess \
    wget git flex bison libgmp3-dev cmake make g++ gcc

#Install IqTree

RUN cd /usr/local/ \
    && wget https://github.com/iqtree/iqtree3/releases/download/v${IQTREE_VERSION}/iqtree-${IQTREE_VERSION}-Linux-intel.tar.gz \
    && tar -xzvf iqtree-${IQTREE_VERSION}-Linux-intel.tar.gz \
    && mv iqtree-${IQTREE_VERSION}-Linux-intel/bin/iqtree3 /usr/local/bin/iqtree \
    && chmod +x /usr/local/bin/iqtree \
    && rm -rf iqtree-${IQTREE_VERSION}-Linux-intel.tar.gz \
    && rm -rf iqtree-${IQTREE_VERSION}-Linux-intel

# INSTALL RAxML
RUN cd /usr/local/ \
    && wget https://github.com/stamatak/standard-RAxML/archive/v${RAXML_VERSION}.tar.gz \
    && tar -xzvf v${RAXML_VERSION}.tar.gz \
    && rm -rf v${RAXML_VERSION}.tar.gz \
    && cd standard-RAxML-${RAXML_VERSION} \
    && sed -i 's/CC = gcc/CC = gcc -std=gnu17/g' Makefile.PTHREADS.gcc \ 
    && make -f Makefile.PTHREADS.gcc \
    && mv raxmlHPC-PTHREADS /usr/local/bin/raxmlHPC \
    && cd /usr/local \
    && rm -rf /usr/local/standard-RAxML-${RAXML_VERSION}

# INSTALL RAxML-NG
RUN cd /usr/local/ \
    && git clone --recursive https://github.com/amkozlov/raxml-ng \
    && cd raxml-ng \
    && git checkout ${RAXMLNG_VERSION} \
    && mkdir build && cd build \
    && cmake .. \
    && make \
    && mv ../bin/raxml-ng /usr/local/bin/ \
    && cd /usr/local/ \
    && rm -rf raxml-ng

#INSTALL FASTTREE
RUN cd /usr/local \
    && cd /usr/local/ \
    && git clone https://github.com/morgannprice/fasttree.git \
    && cd fasttree \
    && git checkout v${FASTTREE_VERSION} \
    && gcc -DNO_SSE -O3 -finline-functions -funroll-loops -Wall -o /usr/local/bin/FastTree FastTree.c -lm \
    && cd .. && rm -rf fasttree 

#INSTALL RAPIDNJ
RUN cd /usr/local \
    && wget -O rapidnj.tar.gz https://github.com/somme89/rapidNJ/archive/refs/tags/latest.tar.gz \
    && tar -zxvf rapidnj.tar.gz \
    && cd rapidNJ-latest/ \
    && make \
    && mv bin/* /usr/local/bin \
    && cd /usr/local \
    && rm -rf rapidNJ-latest/ rapidnj.tar.gz
    
#INSTALL GUBBINS
RUN apt-get update && apt-get install -y gubbins=3.4.3-2

RUN apt-get autoremove -y \
    && apt-get remove -y git flex bison libgmp3-dev make cmake g++  \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

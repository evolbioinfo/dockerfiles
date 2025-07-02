FROM ubuntu:25.04

LABEL maintainer="Julien FUMEY <julien.fumey@pasteur.fr>"

ENV RAXML_VERSION="8.2.13"
ENV IQTREE_VERSION="3.0.0"
ENV RAXMLNG_VERSION="1.2.2"
ENV FASTTREE_VERSION=""
ENV RAPIDNJ=""

RUN apt-get update --fix-missing && apt-get install -y \
    python3 wget git flex bison libgmp3-dev cmake make g++ gcc

#Install IqTree
RUN cd /usr/local/ \
    && wget https://github.com/iqtree/iqtree3/releases/download/v${IQTREE_VERSION}/iqtree-${IQTREE_VERSION}-Linux-intel.tar.gz \
    && tar -xzvf iqtree-${IQTREE_VERSION}-Linux-intel.tar.gz \
    && mv iqtree-${IQTREE_VERSION}-Linux-intel/bin/iqtree3 /usr/local/bin/iqtree \
    && chmod +x /usr/local/bin/iqtree \
    && cd ../../ \
    && rm -rf iqtree-${IQTREE_VERSION}-Linux-intel.tar.gz \
    && rm -rf iqtree-${IQTREE_VERSION}-Linux-intel

# INSTALL RAxML
RUN cd /usr/local/ \
    && wget https://github.com/stamatak/standard-RAxML/archive/v${RAXML_VERSION}.tar.gz \
    && tar -xzvf v${RAXML_VERSION}.tar.gz \
    && rm -rf v${RAXML_VERSION}.tar.gz \
    && cd standard-RAxML-${RAXML_VERSION} \
    && make -f Makefile.PTHREADS.gcc \
    && mv raxmlHPC-PTHREADS /usr/local/bin/raxmlHPC \
    && cd /usr/local \
    && rm -rf /usr/local/standard-RAxML-${RAXML_VERSION}

# INSTALL RAxML-NG
RUN  apt-get update \
    && cd /usr/local/ \
    && git clone --recursive https://github.com/amkozlov/raxml-ng \
    && cd raxml-ng \
    && git checkout ${RAXMLNG_VERSION} \
    && mkdir build && cd build \
    && cmake .. \
    && make \
    && mv ../bin/raxml-ng /usr/local/bin/ \
    && cd /usr/local/ \
    && rm -rf raxml-ng \
    && apt-get remove -y git flex bison libgmp3-dev make cmake g++ \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

#INSTALL FASTTREE
RUN cd /usr/local \
    && wget https://morgannprice.github.io/fasttree/FastTree.c \
    && gcc -DOPENMP -O3 -fopenmp -fopenmp-simd -funsafe-math-optimizations -march=native -o FastTree FastTree.c -lm \
    && mv FastTree /usr/local/bin \
    && rm FastTree.c

#INSTALL RAPIDNJ
RUN cd /usr/local \
    && wget -o rapidnj.tar.gz https://github.com/somme89/rapidNJ/archive/refs/tags/latest.tar.gz \
    && tar -zxvf rapidnj.tar.gz \
    && cd r
#INSTALL GUBBINS
RUN apt-get update && apt-get install -y gubbins

RUN apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur
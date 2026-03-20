# Gblocks
# Selection of conserved blocks from multiple alignments for their use in phylogenetic analysis 
# Version 1.0
# https://molevol-ibe.csic.es/Gblocks/Gblocks_1.0_Linux.tar.gz

# base image: Ubuntu
FROM ubuntu:24.04

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

RUN apt-get update --fix-missing \
    && apt-get install -y wget \ 
    && cd /usr/local/ \
    && wget https://molevol-ibe.csic.es/Gblocks/Gblocks_1.0_Linux.tar.gz \
    && gunzip Gblocks_1.0_Linux.tar.gz \
    && tar -xvf Gblocks_1.0_Linux.tar \
    && rm -rf Gblocks_1.0_Linux.tar.gz  \
    && cd Gblocks_1.0_Linux \
    && mv Gblocks /usr/local/bin/Gblocks \
    && cd /usr/local \
    && rm -rf Gblocks_1.0_Linux \
    && apt-get remove -y wget \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/Gblocks"]

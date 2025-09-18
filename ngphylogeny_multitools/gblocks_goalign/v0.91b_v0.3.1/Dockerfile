# Gblocks
# Selection of conserved blocks from multiple alignments for their use in phylogenetic analysis 
# http://molevol.cmima.csic.es/castresana/Gblocks.html
# Version 0.91b, January 2002
# http://molevol.cmima.csic.es/castresana/Gblocks/Gblocks_Linux64_0.91b.tar.Z
# goalign
# https://github.com/evolbioinfo/goalign
# v0.3.1

# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Frédéric Lemoine <frederic.lemoine@pasteur.fr>

ENV PATH=/usr/local/go/bin:/gopath/bin/:$PATH
ENV GOPATH=/gopath

RUN apt-get update --fix-missing \
    && apt-get install -y wget \ 
    && cd /usr/local/ \
    && wget --no-check-certificate -O Gblocks_Linux64_0.91b.tar.Z https://molevol-ibe.csic.es/Gblocks/Gblocks_Linux64_0.91b.tar.Z \
    && tar -xzvf Gblocks_Linux64_0.91b.tar.Z \
    && rm -rf Gblocks_Linux64_0.91b.tar.Z \
    && cd Gblocks_0.91b \
    && mv Gblocks /usr/local/bin/Gblocks \
    && cd /usr/local \
    && rm -rf Gblocks_0.91b

RUN cd /usr/local/bin \
    && wget https://github.com/evolbioinfo/goalign/releases/download/v0.3.1/goalign_amd64_linux \
    && mv goalign_amd64_linux goalign \
    && chmod +x goalign \
    && apt-get remove -y wget \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/Gblocks"]

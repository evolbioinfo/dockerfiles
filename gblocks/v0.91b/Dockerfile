# Gblocks
# Selection of conserved blocks from multiple alignments for their use in phylogenetic analysis 
# http://molevol.cmima.csic.es/castresana/Gblocks.html
# Version 0.91b, January 2002
# http://molevol.cmima.csic.es/castresana/Gblocks/Gblocks_Linux64_0.91b.tar.Z

# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Mathieu Moslonka-Lefebvre <mathieu.moslonka-lefebvre@pasteur.fr>

RUN apt-get update --fix-missing \
    && apt-get install -y wget \ 
    && cd /usr/local/ \
    && wget -O Gblocks_Linux64_0.91b.tar.Z http://molevol.cmima.csic.es/castresana/Gblocks/Gblocks_Linux64_0.91b.tar.Z \
    && tar -xzvf Gblocks_Linux64_0.91b.tar.Z \
    && rm -rf Gblocks_Linux64_0.91b.tar.Z \
    && cd Gblocks_0.91b \
    && mv Gblocks /usr/local/bin/Gblocks \
    && cd /usr/local \
    && rm -rf Gblocks_0.91b \
    && apt-get remove -y wget \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/Gblocks"]

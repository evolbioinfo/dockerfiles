# IQ-TREE
# http://www.iqtree.org/
# 2.0.6

# base image: Ubuntu
FROM ubuntu:20.04

# File Author / Maintainer
LABEL maintainer=frederic.lemoine@pasteur.fr

RUN apt-get update --fix-missing \
    && apt-get install -y wget \
    && cd /usr/local/ \
    && wget https://github.com/Cibiv/IQ-TREE/releases/download/v2.0.6/iqtree-2.0.6-Linux.tar.gz \
    && tar -xzvf iqtree-2.0.6-Linux.tar.gz \
    && mv iqtree-2.0.6-Linux/bin/iqtree2 /usr/local/bin/iqtree \
    && chmod +x /usr/local/bin/iqtree \
    && cd ../../ \
    && rm -rf iqtree-2.0.6-Linux.tar.gz \
    && rm -rf iqtree-2.0.6-Linux \
    && apt-get remove -y wget \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["iqtree"]

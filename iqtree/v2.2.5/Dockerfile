# IQ-TREE
# http://www.iqtree.org/
# 2.2.5

# base image: Ubuntu
FROM ubuntu:22.04

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

ENV VERSION="2.2.5"

RUN apt-get update --fix-missing \
    && apt-get install -y wget \
    && cd /usr/local/ \
    && wget https://github.com/iqtree/iqtree2/releases/download/v${VERSION}/iqtree-${VERSION}-Linux.tar.gz \
    && tar -xzvf iqtree-${VERSION}-Linux.tar.gz \
    && mv iqtree-${VERSION}-Linux/bin/iqtree2 /usr/local/bin/iqtree \
    && chmod +x /usr/local/bin/iqtree \
    && cd ../../ \
    && rm -rf iqtree-${VERSION}-Linux.tar.gz \
    && rm -rf iqtree-${VERSION}-Linux \
    && apt-get remove -y wget \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["iqtree"]

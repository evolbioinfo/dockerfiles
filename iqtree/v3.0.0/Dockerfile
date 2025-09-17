# IQ-TREE
# http://www.iqtree.org/
# 3.0.0

# base image: Ubuntu
FROM ubuntu:24.04

# File Author / Maintainer
LABEL maintainer="julien.fumey@pasteur.fr"

ENV VERSION="3.0.0"

RUN apt-get update --fix-missing \
    && apt-get install -y wget \
    && cd /usr/local/ \
    && wget https://github.com/iqtree/iqtree3/releases/download/v${VERSION}/iqtree-${VERSION}-Linux-intel.tar.gz \
    && tar -xzvf iqtree-${VERSION}-Linux-intel.tar.gz \
    && mv iqtree-${VERSION}-Linux-intel/bin/iqtree3 /usr/local/bin/iqtree \
    && chmod +x /usr/local/bin/iqtree \
    && cd ../../ \
    && rm -rf iqtree-${VERSION}-Linux-intel.tar.gz \
    && rm -rf iqtree-${VERSION}-Linux-intel \
    && apt-get remove -y wget \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["iqtree"]
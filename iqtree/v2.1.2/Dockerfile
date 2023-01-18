# IQ-TREE
# http://www.iqtree.org/
# 2.1.2

# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Anna Zhukova <anna.zhukova@pasteur.fr>

RUN apt-get update --fix-missing \
    && apt-get install -y wget \
    && cd /usr/local/ \
    && wget https://github.com/iqtree/iqtree2/releases/download/v2.1.2/iqtree-2.1.2-Linux.tar.gz \
    && tar -xzvf iqtree-2.1.2-Linux.tar.gz \
    && mv iqtree-2.1.2-Linux/bin/iqtree2 /usr/local/bin/iqtree \
    && chmod +x /usr/local/bin/iqtree \
    && cd ../../ \
    && rm -rf iqtree-2.1.2-Linux.tar.gz \
    && rm -rf iqtree-2.1.2-Linux \
    && apt-get remove -y wget \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["iqtree"]

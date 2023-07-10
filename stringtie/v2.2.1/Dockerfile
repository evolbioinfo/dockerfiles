# stringtie
# https://ccb.jhu.edu/software/stringtie/
# v2.2.1

# base image: Ubuntu
FROM ubuntu:22.04

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

ENV VERSION=2.2.1


RUN apt-get update --fix-missing \
    && apt-get install -y wget python2 \
    && wget  http://ccb.jhu.edu/software/stringtie/dl/stringtie-2.2.1.Linux_x86_64.tar.gz \
    && tar -xzvf stringtie-2.2.1.Linux_x86_64.tar.gz \
    && mv stringtie-2.2.1.Linux_x86_64/* /usr/local/bin/ \
    && rm -rf stringtie-2.2.1.Linux_x86_64 \
    && apt-get remove -y wget gcc make libbz2-dev zlib1g-dev liblzma-dev libncurses-dev bzip2 \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/stringtie"]

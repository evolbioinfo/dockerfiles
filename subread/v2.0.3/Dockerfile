# Subread
# https://sourceforge.net/projects/subread/
# 2.0.3

# base image: Ubuntu
FROM ubuntu:22.04

ENV VERSION=2.0.3


# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

RUN apt-get update --fix-missing \
    && apt-get install -y wget \
    && cd /usr/local/ \
    && wget -O subread-${VERSION}-Linux-x86_64.tar.gz "https://downloads.sourceforge.net/project/subread/subread-${VERSION}/subread-${VERSION}-Linux-x86_64.tar.gz?r=https%3A%2F%2Fsourceforge.net%2Fprojects%2Fsubread%2Ffiles%2Fsubread-${VERSION}%2Fsubread-${VERSION}-Linux-x86_64.tar.gz%2Fdownload%3Fuse_mirror%3Dfreefr&ts=1601456911&use_mirror=freefr" \
    && tar -xzvf subread-${VERSION}-Linux-x86_64.tar.gz \
    && mv subread-${VERSION}-Linux-x86_64/bin/* /usr/local/bin/ \
    && mv subread-${VERSION}-Linux-x86_64/annotation /usr/local/ \
    && rm -rf subread-${VERSION}-Linux-x86_64 subread-${VERSION}-Linux-x86_64.tar.gz \
    && apt-get remove -y wget \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENV PATH=/usr/local/bin/utilities:$PATH

ENTRYPOINT ["/usr/local/bin/featureCounts"]

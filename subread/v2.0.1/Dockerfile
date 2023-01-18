# Subread
# https://sourceforge.net/projects/subread/
# 2.0.1

# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

RUN apt-get update --fix-missing \
    && apt-get install -y wget \
    && cd /usr/local/ \
    && wget -O subread-2.0.1-Linux-x86_64.tar.gz "https://downloads.sourceforge.net/project/subread/subread-2.0.1/subread-2.0.1-Linux-x86_64.tar.gz?r=https%3A%2F%2Fsourceforge.net%2Fprojects%2Fsubread%2Ffiles%2Fsubread-2.0.1%2Fsubread-2.0.1-Linux-x86_64.tar.gz%2Fdownload%3Fuse_mirror%3Dfreefr&ts=1601456911&use_mirror=freefr" \
    && tar -xzvf subread-2.0.1-Linux-x86_64.tar.gz \
    && mv subread-2.0.1-Linux-x86_64/bin/* /usr/local/bin/ \
    && mv subread-2.0.1-Linux-x86_64/annotation /usr/local/ \
    && rm -rf subread-2.0.1-Linux-x86_64 subread-2.0.1-Linux-x86_64.tar.gz \
    && apt-get remove -y wget \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENV PATH=/usr/local/bin/utilities:$PATH

ENTRYPOINT ["/usr/local/bin/featureCounts"]

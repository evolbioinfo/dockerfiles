# Seq-Gen
# https://github.com/rambaut/Seq-Gen/
# http://tree.bio.ed.ac.uk/software/seqgen/
# v1.3.4_HIVB

# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

RUN apt-get update --fix-missing \
    && apt-get install -y wget make gcc \
    && cd /usr/local/ \
    && wget -O Seq-Gen-1.3.4_hivb.tar.gz https://github.com/fredericlemoine/Seq-Gen/archive/1.3.4_hivb.tar.gz \
    && tar xvzf Seq-Gen-1.3.4_hivb.tar.gz \
    && rm -f Seq-Gen-1.3.4_hivb.tar.gz \
    && cd /usr/local/Seq-Gen-1.3.4_hivb/source \
    && make \
    && mv seq-gen /usr/local/bin/ \
    && cd /usr/local \
    && rm -rf Seq-Gen-1.3.4_hivb \
    && apt-get remove -y wget make gcc \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/seq-gen"]

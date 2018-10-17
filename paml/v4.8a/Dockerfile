# PAML
# http://abacus.gene.ucl.ac.uk/software/paml.html
# v4.8a

# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

RUN apt-get update --fix-missing \
    && apt-get install -y wget gcc make \
    && cd /usr/local/ \
    && wget http://abacus.gene.ucl.ac.uk/software/paml4.8a.macosx.tgz \
    && tar -xzvf paml4.8a.macosx.tgz \
    && rm -rf paml4.8a.macosx.tgz \
    && cd paml4.8/src \
    && make \
    && cp baseml basemlg chi2 codeml evolver infinitesites mcmctree pamp yn00 /usr/local/bin/ \
    && cd /usr/local \
    && rm -rf /usr/local/paml4.8 \
    && apt-get remove -y wget gcc make \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["codeml"]

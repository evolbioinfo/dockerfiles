# bowtie
# https://github.com/BenLangmead/bowtie/
# v0.12.9

# base image: Ubuntu
FROM ubuntu:22.04

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

ENV VERSION=v0.12.9

RUN apt-get update --fix-missing \
    && apt-get install -y git g++ make \
    && cd /usr/local/ \
    && git clone https://github.com/BenLangmead/bowtie.git \
    && cd bowtie \
    && git checkout $VERSION \
    && make \
    && cp bowtie bowtie-build bowtie-inspect /usr/local/bin/ \
    && cd /usr/local \
    && rm -rf bowtie \
    && apt-get remove -y git gcc make \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/bowtie"]

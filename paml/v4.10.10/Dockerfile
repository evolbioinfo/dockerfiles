# PAML
# https://github.com/abacus-gene/paml
# v4.10.10

# base image: Ubuntu
FROM ubuntu:24.04

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

ENV VERSION=4.10.10

# Linux binaries published by upstream; data files (models, matrices)
# are kept in /usr/local/share/paml
RUN apt-get update --fix-missing \
    && apt-get install -y wget \
    && cd /usr/local/ \
    && wget https://github.com/abacus-gene/paml/releases/download/v${VERSION}/paml-${VERSION}-linux-x86_64.tar.gz \
    && tar -xzvf paml-${VERSION}-linux-x86_64.tar.gz \
    && rm -f paml-${VERSION}-linux-x86_64.tar.gz \
    && chmod +x paml-${VERSION}-linux-x86_64/bin/* \
    && mv paml-${VERSION}-linux-x86_64/bin/* /usr/local/bin/ \
    && mkdir -p /usr/local/share/paml \
    && mv paml-${VERSION}-linux-x86_64/dat /usr/local/share/paml/ \
    && rm -rf paml-${VERSION}-linux-x86_64 \
    && apt-get remove -y wget \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["codeml"]

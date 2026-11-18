# adapter removal
# https://github.com/MikkelSchubert/adapterremoval
# v3.0.2
# Note: AdapterRemoval 3 renames the executable to adapterremoval3 and changes
# its command line, see https://adapterremoval.readthedocs.io/en/latest/migrating.html

# base image: Ubuntu
FROM ubuntu:24.04

# File Author / Maintainer
LABEL maintainer=frederic.lemoine@pasteur.fr

ENV VERSION=3.0.2

# Statically linked binary published by upstream
RUN apt-get update --fix-missing \
    && apt-get install -y wget \
    && cd /usr/local/ \
    && wget https://github.com/MikkelSchubert/adapterremoval/releases/download/v${VERSION}/adapterremoval-v${VERSION}-linux-x86_64.tar.gz \
    && tar -xzvf adapterremoval-v${VERSION}-linux-x86_64.tar.gz \
    && mv adapterremoval-v${VERSION}-linux-x86_64/adapterremoval3 /usr/local/bin/ \
    && chmod +x /usr/local/bin/adapterremoval3 \
    && rm -rf adapterremoval-v${VERSION}-linux-x86_64 adapterremoval-v${VERSION}-linux-x86_64.tar.gz \
    && apt-get remove -y wget \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/adapterremoval3"]

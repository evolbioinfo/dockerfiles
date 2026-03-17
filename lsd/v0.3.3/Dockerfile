# LSD
# https://github.com/tothuhien/lsd-0.3beta
# v0.3.3 - revision 55183ca9d0ea918314e7c89e3dfd67382b1f2597
# https://github.com/tothuhien/lsd-0.3beta/

# base image: Ubuntu
FROM ubuntu:24.04

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

RUN apt-get update --fix-missing \
    && apt-get install -y wget \ 
    && cd /usr/local/ \
    && wget https://github.com/tothuhien/lsd-0.3beta/releases/download/v0.3.3/lsd_unix \
    && mv lsd_unix /usr/local/bin/lsd \
    && chmod +x /usr/local/bin/lsd \
    && apt-get remove -y wget \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/lsd"]

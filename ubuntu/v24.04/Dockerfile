# Ubuntu
# https://www.ubuntu.com/
# v24.04

# base image: ubuntu:24.04
FROM ubuntu:24.04

# File Author / Maintainer
LABEL maintainer=frederic.lemoine@pasteur.fr


RUN apt-get update --fix-missing \
    && apt-get install -y curl wget pigz xz-utils perl \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/bin/bash"]

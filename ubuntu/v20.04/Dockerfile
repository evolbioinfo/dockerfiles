# Ubuntu
# https://www.ubuntu.com/
# v20.04

# base image: ubuntu:20.04
FROM ubuntu:20.04

# File Author / Maintainer
LABEL maintainer=frederic.lemoine@pasteur.fr


RUN apt-get update --fix-missing \
    && apt-get install -y curl wget pigz \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/bin/bash"]

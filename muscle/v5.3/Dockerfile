# MUSCLE
# http://www.drive5.com/muscle/downloads.htm
# v5.3

# base image: Ubuntu
FROM ubuntu:24.04

ENV VERSION="5.3"

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

RUN apt-get update --fix-missing \
    && apt-get install -y wget git make g++ file python3 curl \
    && cd /usr/local/ \
    && git clone https://github.com/rcedgar/muscle.git \
    && cd muscle \
    && git checkout v$VERSION \
    && cd src \
    && curl -fsSL https://raw.githubusercontent.com/rcedgar/vcxproj_make/806d016/vcxproj_make.py > vcxproj_make.py \
    && chmod +x vcxproj_make.py \
    && ./vcxproj_make.py --openmp \
    && mv ../bin/muscle /usr/local/bin \
    && cd /usr/local \
    && rm -rf muscle \
    && apt-get remove -y wget g++ make curl \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["muscle"]

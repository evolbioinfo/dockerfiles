# BMGE
# https://gitlab.pasteur.fr/GIPhy/BMGE
# v2.00

# base image: ubuntu:24.04
FROM ubuntu:24.04

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

ENV VERSION="2.00"

COPY BMGE.sh /usr/local/bin/BMGE

RUN apt-get update --fix-missing \
    && apt-get install -y wget git openjdk-25-jre  openjdk-25-jdk \
    && cd /usr/local/ \
    && git clone https://gitlab.pasteur.fr/GIPhy/BMGE.git  \
    && cd BMGE \
    && git checkout ${VERSION} \
    && cd src \
    && javac BMGE.java \
    && echo Main-Class: BMGE > MANIFEST.MF \
    && jar -cmvf MANIFEST.MF BMGE.jar BMGE.class bmge/*.class \
    && mv BMGE.jar /usr/local/bin/ \
    && chmod +x /usr/local/bin/BMGE \
    && apt-get remove -y wget openjdk-25-jdk git \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/BMGE"]

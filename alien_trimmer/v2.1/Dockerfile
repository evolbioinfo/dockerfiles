# Alien Trimmer
# ftp://ftp.pasteur.fr/pub/gensoft/projects/AlienTrimmer/
# v2.1

# base image: ubuntu:16.04
FROM ubuntu:22.04

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

COPY AlienTrimmer.sh /usr/local/bin/AlienTrimmer

RUN apt-get update --fix-missing \
    && apt-get install -y wget git openjdk-11-jre openjdk-11-jdk \
    && cd /usr/local/ \
    && git clone https://gitlab.pasteur.fr/GIPhy/AlienTrimmer.git \
    && cd AlienTrimmer \
    && git checkout 2.1 \
    && cd src \
    && javac AlienTrimmer.java \
    && chmod +x /usr/local/bin/AlienTrimmer \
    && apt-get remove -y wget openjdk-11-jdk \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/AlienTrimmer"]

# Alien Trimmer
# ftp://ftp.pasteur.fr/pub/gensoft/projects/AlienTrimmer/
# v0.4.0

# base image: ubuntu:16.04
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

COPY AlienTrimmer.sh /usr/local/bin/AlienTrimmer

RUN apt-get update --fix-missing \
    && apt-get install -y wget openjdk-8-jre openjdk-8-jdk \
    && cd /usr/local/ \
    && wget ftp://ftp.pasteur.fr/pub/gensoft/projects/AlienTrimmer/AlienTrimmer_0.4.0.tar.gz \
    && tar -xzvf AlienTrimmer_0.4.0.tar.gz \
    && rm -rf AlienTrimmer_0.4.0.tar.gz \
    && cd AlienTrimmer_0.4.0/src \
    && javac AlienTrimmer.java \
    && chmod +x /usr/local/bin/AlienTrimmer \
    && apt-get remove -y wget openjdk-8-jdk \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/AlienTrimmer"]

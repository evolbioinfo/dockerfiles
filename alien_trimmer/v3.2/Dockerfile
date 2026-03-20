# Alien Trimmer
# https://gitlab.pasteur.fr/GIPhy/AlienTrimmer
# v3.2

# base image: ubuntu:16.04
FROM ubuntu:22.04

ENV VERSION="3.2"

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

RUN apt-get update --fix-missing \
    && apt-get install -y wget git openjdk-25-jre openjdk-25-jdk \
    && cd /usr/local/ \
    && git clone https://gitlab.pasteur.fr/GIPhy/AlienTrimmer.git \
    && cd AlienTrimmer \
    && git checkout ${VERSION} \
    && cd src \
    && javac AlienTrimmer.java \
    && echo Main-Class: AlienTrimmer > MANIFEST.MF \
    && jar -cmvf MANIFEST.MF AlienTrimmer.jar AlienTrimmer*.class \
    && rm MANIFEST.MF AlienTrimmer*.class \
    && mv AlienTrimmer.sh /usr/local/bin/AlienTrimmer \
    && mv AlienTrimmer.jar /usr/local/bin/ \
    && chmod +x /usr/local/bin/AlienTrimmer \
    && apt-get remove -y wget openjdk-25-jdk \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/AlienTrimmer"]

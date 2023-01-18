# haplogrep
# https://github.com/seppinho/haplogrep-cmd
# v2.4.0

# base image: Ubuntu
FROM ubuntu:20.04

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

ENV VERSION=v2.4.0

ENV PATH $PATH:/usr/local/IGV_${FULLVERSION}
ENV DEBIAN_FRONTEND="noninteractive" TZ="Europe/Paris"

RUN apt-get update --fix-missing \
    && apt-get install -y unzip default-jdk xvfb wget \
    && cd /usr/local/bin \
    && wget https://github.com/seppinho/haplogrep-cmd/releases/download/${VERSION}/haplogrep.zip \
    && unzip haplogrep.zip \
    && rm haplogrep.zip \
    && chmod +x haplogrep \
    && apt-get remove -y unzip \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/haplogrep"]

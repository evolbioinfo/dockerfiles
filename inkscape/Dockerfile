# inkscape
# https://inkscape.org/fr/
# Stable

# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

RUN apt-get update --fix-missing \
    && echo "deb http://ppa.launchpad.net/inkscape.dev/stable/ubuntu xenial main" >> /etc/apt/sources.list \
    && apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 22C4850146603F3DB0ED00479DA4BD18B9A06DE3 \
    && apt-get update \
    && apt-get install -y inkscape \
    && apt-get install ttf-ubuntu-font-family \
    && mkdir /pasteur

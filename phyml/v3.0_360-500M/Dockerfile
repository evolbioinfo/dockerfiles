# PhyML
# http://www.atgc-montpellier.fr/phyml/
# Downloaded from https://www.softpedia.com/get/Science-CAD/PhyML.shtml
# v3.0-360-500M

# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

RUN apt-get update --fix-missing \
    && apt-get install -y wget \
    && wget http://www.atgc-montpellier.fr/download/binaries/phyml/phyml_v3.0.tar.gz \
    && tar -xzvf phyml_v3.0.tar.gz \
    && mv PhyML_3.0/PhyML_3.0_linux64 /usr/local/bin/phyml \
    && rm -rf PhyML_3.0/ phyml_v3.0.tar.gz \
    && apt-get remove -y wget \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["phyml"]

# IGV
# http://software.broadinstitute.org/software/igv/
# https://github.com/igvteam/igv
# v2.19.8

# base image: Ubuntu
FROM ubuntu:24.04

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

ENV VERSION=2.19
ENV FULLVERSION=2.19.8
ENV PATH=$PATH:/usr/local/IGV_${FULLVERSION}
ENV DEBIAN_FRONTEND="noninteractive" TZ="Europe/Paris"

# IGV 2.19 requires Java 21
RUN apt-get update --fix-missing \
    && apt-get install -y unzip openjdk-21-jre xvfb wget xorg \
    && cd /usr/local/ \
    && wget https://data.broadinstitute.org/igv/projects/downloads/${VERSION}/IGV_${FULLVERSION}.zip \
    && unzip IGV_${FULLVERSION}.zip \
    && rm -f IGV_${FULLVERSION}.zip \
    && cd IGV_${FULLVERSION} \
    && mv igv.sh igv \
    && apt-get remove -y unzip wget \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/IGV_2.19.8/igv"]

# IGV
# http://software.broadinstitute.org/software/igv/
# v2.9

# base image: Ubuntu
FROM ubuntu:20.04

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

ENV VERSION=2.9
ENV FULLVERSION=2.9.0


ENV PATH $PATH:/usr/local/IGV_${FULLVERSION}
ENV DEBIAN_FRONTEND="noninteractive" TZ="Europe/Paris"

RUN apt-get update --fix-missing \
    && apt-get install -y unzip default-jdk xvfb wget xorg \
    && cd /usr/local/ \
    && wget https://data.broadinstitute.org/igv/projects/downloads/${VERSION}/IGV_${FULLVERSION}.zip \
    && unzip IGV_${FULLVERSION}.zip \
    && cd IGV_${FULLVERSION} \
    && mv igv.sh igv \
    && apt-get remove -y unzip \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/IGV_Linux_$FULLVERSION/igv"]

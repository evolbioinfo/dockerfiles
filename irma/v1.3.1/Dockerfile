# IRMA
# https://wonder.cdc.gov/amd/flu/irma/
# https://github.com/CDCgov/irma
# v1.3.1

# base image: r
FROM ubuntu:24.04

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

ENV DEBIAN_FRONTEND="noninteractive" TZ="Europe/Paris"

RUN apt-get update --fix-missing \
    && apt-get install -y zip perl wget r-base \
    && wget https://github.com/CDCgov/irma/releases/download/v1.3.1/irma-v1.3.1-universal.zip \
    && unzip irma-v1.3.1-universal.zip \
    && cd flu-amd/ \
    && mv * /usr/local/bin/ \
    && cd ../ \
    && rm -rf irma-v1.3.1-universal.zip flu-amd \
    && cat /usr/local/bin/IRMA_RES/modules/FLU/config/FLU-minion.sh > /usr/local/bin/IRMA_RES/modules/FLU/config/FLU-miniontmp.sh \
    && echo "TMP=\$TMPDIR" >> /usr/local/bin/IRMA_RES/modules/FLU/config/FLU-miniontmp.sh \
    && cat /usr/local/bin/IRMA_RES/modules/FLU/config/FLU.sh > /usr/local/bin/IRMA_RES/modules/FLU/config/FLU-tmp.sh \
    && sed -i "$ a TMP=\$TMPDIR"  /usr/local/bin/IRMA_RES/modules/*/config/*.sh\
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/IRMA"]

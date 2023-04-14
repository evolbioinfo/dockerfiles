# IRMA
# https://wonder.cdc.gov/amd/flu/irma/
# v1.0.3

# base image: r
FROM ubuntu:22.04

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

ENV DEBIAN_FRONTEND="noninteractive" TZ="Europe/Paris"

RUN apt-get update --fix-missing \
    && apt-get install -y zip perl wget r-base \
    && wget https://wonder.cdc.gov/amd/flu/irma/flu-amd-202209.zip \
    && unzip flu-amd-202209.zip \
    && cd flu-amd/ \
    && mv * /usr/local/bin/ \
    && cd ../ \
    && rm -rf flu-amd-201909.zip flu-amd \
    && cat /usr/local/bin/IRMA_RES/modules/FLU/config/FLU-minion.sh > /usr/local/bin/IRMA_RES/modules/FLU/config/FLU-miniontmp.sh \
    && echo "TMP=\$TMPDIR" >> /usr/local/bin/IRMA_RES/modules/FLU/config/FLU-miniontmp.sh \
    && cat /usr/local/bin/IRMA_RES/modules/FLU/config/FLU.sh > /usr/local/bin/IRMA_RES/modules/FLU/config/FLU-tmp.sh \
    && sed -i "$ a TMP=\$TMPDIR"  /usr/local/bin/IRMA_RES/modules/*/config/*.sh\
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/IRMA"]

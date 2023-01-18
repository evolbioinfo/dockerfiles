# IRMA
# https://wonder.cdc.gov/amd/flu/irma/
# v0.6.7

# base image: r
FROM ubuntu:xenial

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

RUN apt-get update --fix-missing \
    && apt-get install -y zip perl wget r-base \
    && wget https://wonder.cdc.gov/amd/flu/irma/flu-amd-201704.zip \
    && unzip flu-amd-201704.zip \
    && cd flu-amd/ \
    && mv * /usr/local/bin/ \
    && cd ../ \
    && rm -rf flu-amd-201704.zip flu-amd \
    && cat /usr/local/bin/IRMA_RES/modules/FLU/config/FLU-minion.sh > /usr/local/bin/IRMA_RES/modules/FLU/config/FLU-miniontmp.sh \
    && echo "TMP=\$TMPDIR" >> /usr/local/bin/IRMA_RES/modules/FLU/config/FLU-miniontmp.sh \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/IRMA"]

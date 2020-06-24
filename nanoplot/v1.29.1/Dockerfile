# NanoPlot
# https://github.com/wdecoster/NanoPlot
# v1.29.1

# base image: Ubuntu
FROM evolbioinfo/python:v3.8.2

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

ENV VERSION=1.29.1

# Switch to your new user in the docker image
USER root

# Install treetime
RUN cd /usr/local/ && pip3 install --no-cache-dir NanoPlot==${VERSION}

# Switch to your new user in the docker image
USER evolbioinfo

# The entrypoint runs pyhton3
ENTRYPOINT ["NanoPlot"]

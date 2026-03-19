FROM evolbioinfo/python:v3.6.9

# Switch to your new user in the docker image
USER root

ENV VERSION="0.2.27"

# Install bc
RUN apt-get update -y && apt-get install -y bc

# Install treesimulator
RUN cd /usr/local/ && pip3 install --no-cache-dir treesimulator==${VERSION}

# Switch to your new user in the docker image
USER evolbioinfo

# start working in the "evolbioinfo" home directory
WORKDIR /home/evolbioinfo

# The entrypoint runs command line with command line arguments
ENTRYPOINT ["/bin/bash"]

FROM evolbioinfo/python:v3.6.9

# Switch to your new user in the docker image
USER root

# Copy in the version files
ADD *_VERSION /usr/local/

# Install treetime
RUN cd /usr/local/ && pip3 install --no-cache-dir pastml==`cat PASTML_VERSION`

# Switch to your new user in the docker image
USER evolbioinfo

# start working in the "evolbioinfo" home directory
WORKDIR /home/evolbioinfo

# The entrypoint runs pastml with command line arguments
ENTRYPOINT ["pastml"]

FROM evolbioinfo/python:v3.8.2

# Switch to your new user in the docker image
USER root

# Install pastml
RUN cd /usr/local/ && pip3 install --no-cache-dir pastml==1.9.40

# Switch to your new user in the docker image
USER evolbioinfo

# start working in the "evolbioinfo" home directory
WORKDIR /home/evolbioinfo

# The entrypoint runs pastml with command line arguments
ENTRYPOINT ["pastml"]

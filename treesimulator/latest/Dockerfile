FROM evolbioinfo/python:v3.6.9

# Switch to your new user in the docker image
USER root

# Install treesimulator
RUN cd /usr/local/ && pip3 install --no-cache-dir treesimulator==0.1.2

# Switch to your new user in the docker image
USER evolbioinfo

# start working in the "evolbioinfo" home directory
WORKDIR /home/evolbioinfo

# The entrypoint runs command line with command line arguments
ENTRYPOINT ["/bin/bash"]
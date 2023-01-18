FROM evolbioinfo/python:v3.6.9

# Switch to your new user in the docker image
USER root

# Install jphmm_tools
RUN cd /usr/local/ && pip3 install --no-cache-dir jphmm_tools==0.1.4

# Switch to your new user in the docker image
USER evolbioinfo

# start working in the "evolbioinfo" home directory
WORKDIR /home/evolbioinfo

# The entrypoint runs bash
ENTRYPOINT ["/bin/bash"]

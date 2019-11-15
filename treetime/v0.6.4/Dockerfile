# TreeTime
# https://github.com/neherlab/treetime
# v0.6.4

FROM evolbioinfo/python:v3.6.9

# Switch to your new user in the docker image
USER root

# Install treetime
RUN cd /usr/local/ && pip3 install --no-cache-dir phylo-treetime==0.6.4

# Switch to your new user in the docker image
USER evolbioinfo

# The entrypoint runs treetime with command line arguments
ENTRYPOINT ["treetime"]


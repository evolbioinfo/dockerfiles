# Version 0.9

#FROM condaforge/miniforge3:26.1.0-0
FROM continuumio/miniconda3:4.12.0

ENV VERSION=0.9

# Install phylodeep
RUN cd /usr/local/ && pip3 install --no-cache-dir phylodeep==${VERSION}

# The entrypoint runs command line with command line argument
ENTRYPOINT ["/bin/bash"]
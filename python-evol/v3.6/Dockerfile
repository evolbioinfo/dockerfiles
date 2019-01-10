# python3.6 with packages:
# - ete3==3.1.1
# - biopython==1.72
# - numpy==1.15.4
# - pandas==0.23.4
# - scipy==1.1.0
# - matplotlib==3.0.2
# - hdx-python-country==1.5.4

FROM python:3.6-slim

# Python3 needs an UTF-8 locale, http://bugs.python.org/issue19846
ENV LANG C.UTF-8

# Build-time environmental variable so that apt doesn't complain
ARG DEBIAN_FRONTEND=noninteractive

# File Author / Maintainer
MAINTAINER Anna Zhukova <anna.zhukova@pasteur.fr>

# Install pastml
RUN cd /usr/local/ && pip3 install --no-cache-dir ete3==3.1.1 biopython==1.72 pandas==0.23.4 scipy==1.1.0 matplotlib==3.0.2 hdx-python-country==1.5.4 numpy==1.15.4

# Clean up
RUN mkdir /pasteur

# create a non-root user named evolbioinfo,
# give them the password "evolbioinfo" put them in the sudo group
RUN useradd -d /home/evolbioinfo -m -s /bin/bash evolbioinfo && echo "evolbioinfo:evolbioinfo" | chpasswd && adduser evolbioinfo sudo

# start working in the "evolbioinfo" home directory
WORKDIR /home/evolbioinfo

# Make the files owned by evolbioinfo
RUN chown -R evolbioinfo:evolbioinfo /home/evolbioinfo

# Switch to your new user in the docker image
USER evolbioinfo

# The entrypoint runs pyhton3
ENTRYPOINT ["python3"]

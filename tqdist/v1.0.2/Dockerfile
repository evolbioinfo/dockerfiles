FROM evolbioinfo/gotree:v0.4.4

# Python3 needs an UTF-8 locale, http://bugs.python.org/issue19846
ENV LANG C.UTF-8

# Build-time environmental variable so that apt doesn't complain
ARG DEBIAN_FRONTEND=noninteractive
# File Author / Maintainer
LABEL maintainer=frederic.lemoine@pasteur.fr


# Create the development environment
RUN apt clean && apt-get update --fix-missing \
    && apt install -y build-essential cmake wget \
    && cd /usr/local/ \
    && wget https://www.birc.au.dk/~cstorm/software/tqdist/files/tqDist-1.0.2.tar.gz \
    && tar -xvf tqDist-1.0.2.tar.gz \
    && rm -rf tqDist-1.0.2.tar.gz \
    && cd tqDist-1.0.2 \
    && cmake . \
    && make -j 4 \
    && make install \
    && cd /usr/local \
    && rm -rf /usr/local/tqDist-1.0.2 \
    && apt-get remove build-essential -y \
    && apt-get autoremove -y \
    && apt-get autoclean -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir -p /pasteur

# create a non-root user named evolbioinfo,
# give them the password "evolbioinfo" put them in the sudo group
RUN useradd -d /home/evolbioinfo -m -s /bin/bash evolbioinfo && echo "evolbioinfo:evolbioinfo" | chpasswd && adduser evolbioinfo sudo

# start working in the "evolbioinfo" home directory
WORKDIR /home/evolbioinfo

# Make the files owned by evolbioinfo
RUN chown -R evolbioinfo:evolbioinfo /home/evolbioinfo

# Switch to your new user in the docker image
USER evolbioinfo

# The entrypoint runs quartet_dist with command line arguments
ENTRYPOINT ["/usr/local/bin/quartet_dist"]

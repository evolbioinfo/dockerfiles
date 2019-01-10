FROM ubuntu:trusty

# Python3 needs an UTF-8 locale, http://bugs.python.org/issue19846
ENV LANG C.UTF-8

# Build-time environmental variable so that apt doesn't complain
ARG DEBIAN_FRONTEND=noninteractive

# Create the development environment
RUN apt-get update --fix-missing \
    && apt install -y build-essential cmake wget \
    && apt-get upgrade -y

# File Author / Maintainer
MAINTAINER Anna Zhukova <anna.zhukova@pasteur.fr>

# Install tqDist (add include of string.h otherwise g++ does not want to compile it)
RUN cd /usr/local/ \
    && wget http://users-cs.au.dk/cstorm/software/tqdist/files/tqDist-1.0.1.tar.gz \
    && tar -xvf tqDist-1.0.1.tar.gz \
    && rm -rf tqDist-1.0.1.tar.gz \
    && cd tqDist-1.0.1 \
    && sed '1i\
\#include <string.h>' tqDist/pairs_quartet_distance.cpp > tqDist/pairs_quartet_distance.cpp.new \
    && mv tqDist/pairs_quartet_distance.cpp.new tqDist/pairs_quartet_distance.cpp \
    && sed '1i\
\#include <string.h>' tqDist/pairs_triplet_distance.cpp  > tqDist/pairs_triplet_distance.cpp.new \
    && mv tqDist/pairs_triplet_distance.cpp.new tqDist/pairs_triplet_distance.cpp \
    && sed '1i\
\#include <string.h>' tqDist/quartet_dist.cpp  > tqDist/quartet_dist.cpp.new \
    && mv tqDist/quartet_dist.cpp.new tqDist/quartet_dist.cpp \
    && sed '1i\
\#include <string.h>' tqDist/triplet_dist.cpp  > tqDist/triplet_dist.cpp.new \
    && mv tqDist/triplet_dist.cpp.new tqDist/triplet_dist.cpp \
    && cmake . \
    && make -j 4 \
    && make install \
    && cd /usr/local \
    && rm -rf /usr/local/tqDist-1.0.1

# Clean up

# Clean up the development environment to reduce image size
RUN apt-get remove build-essential -y \
    && apt-get autoremove -y \
    && apt-get autoclean -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

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
ENTRYPOINT ["/bin/bash"]

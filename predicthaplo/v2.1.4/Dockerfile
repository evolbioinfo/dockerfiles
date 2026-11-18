# predicthaplo
# https://github.com/cbg-ethz/PredictHaplo
# v2.1.4
# Note: PredictHaplo 2.x takes command line options instead of a configuration
# file, and the executable is now called predicthaplo (PredictHaplo is kept as
# a symbolic link).

# base image: Ubuntu
FROM ubuntu:24.04

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

ENV VERSION=2.1.4

# CXXFLAGS is overridden to drop upstream's -Werror, which breaks the build
# with recent compilers
RUN apt-get update --fix-missing \
    && apt-get install -y wget g++ make libblas-dev liblapack-dev libblas3 liblapack3 \
    && cd /usr/local/ \
    && wget -O predicthaplo.tar.gz https://github.com/cbg-ethz/PredictHaplo/archive/refs/tags/${VERSION}.tar.gz \
    && tar -xzvf predicthaplo.tar.gz \
    && rm -f predicthaplo.tar.gz \
    && cd PredictHaplo-${VERSION} \
    && make CXXFLAGS="-O3 -fPIC -std=c++11" \
    && mv predicthaplo /usr/local/bin/ \
    && ln -s /usr/local/bin/predicthaplo /usr/local/bin/PredictHaplo \
    && cd /usr/local \
    && rm -rf /usr/local/PredictHaplo-${VERSION} \
    && apt-get remove -y wget g++ make libblas-dev liblapack-dev \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/predicthaplo"]

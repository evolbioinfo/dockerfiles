# Schmutzi
# https://github.com/grenaud/schmutzi/
# v1.5.7
# The r-base:4.1.1 base image used up to v1.5.6 can no longer run apt-get
# (expired Debian testing/sid keys) and g++-8 is gone: now based on Ubuntu
# 24.04 with R and the R packages from apt. bamtools v2.5.1 does not build
# with GCC 13 (non-const comparators in Sort.h and
# BamMultiMerger_p.h): patched with sed. v2.5.2
# cannot be used, libgab expects the v2.5.1 header layout.
# Upstream src/Makefile links some tools with ${CXXFLAGS}, which contains -c,
# so they are never linked and 'make install' fails: use $(LDFLAGS) instead.

# base image: Ubuntu
FROM ubuntu:24.04

# File Author / Maintainer
LABEL maintainer=frederic.lemoine@pasteur.fr

ENV SCHMUTZIVERSION=v1.5.7
ENV DEBIAN_FRONTEND="noninteractive" TZ="Europe/Paris"
RUN apt-get update --fix-missing \
    && apt-get install -y --no-install-recommends ca-certificates perl git cmake make g++ zlib1g zlib1g-dev \
        r-base-core r-cran-mass r-cran-fitdistrplus \
    && cd /usr/local/ \
    && git clone https://github.com/grenaud/schmutzi.git \
    && cd schmutzi \
    && git clone https://github.com/pezmaster31/bamtools.git lib/bamtools && cd lib/bamtools && git checkout v2.5.1 \
        && sed -i -E 's/(bool operator\(\)\(const [A-Za-z:]+&[^)]*\))[[:space:]]*$/\1 const/' src/api/algorithms/Sort.h src/api/internal/bam/BamMultiMerger_p.h \
        && cd ../../ \
    && git clone https://github.com/grenaud/libgab.git lib/libgab && cd lib/libgab && git checkout 1.0.4 && cd ../../ \
    && git checkout ${SCHMUTZIVERSION} \
    && sed -i 's/${CXX} ${CXXFLAGS} -o \$@/${CXX} $(LDFLAGS) -o $@/' src/Makefile \
    && make && make install \
    && cd .. \
    && rm -rf schmutzi \
    && apt-get remove -y git cmake make g++ zlib1g-dev \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/bin/schmutzi.pl"]

# treestructure
#
# commit bb0d29f from Apr 11 2025
#
# https://github.com/emvolz-phylodynamics/treestructure
#
# Erik Volz, Carsten Wiuf, Yonatan Grad, Simon Frost, Ann Dennis, Xavier Didelot,
# Systematic Biology, Volume 69, Issue 5, September 2020, Pages 884–896;
# doi: https://doi.org/10.1093/sysbio/syaa009

###############################################################################
# Builder stage
###############################################################################
# base image: r-base:4.6.0
FROM r-base:4.6.0 AS builder

ENV CRANREPO=https://cloud.r-project.org
ENV R_LIBS_SITE=/usr/local/lib/R/site-library

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential     \
    procps              \
    curl                \
    git                 \
    ca-certificates     \
    pandoc              \
    libssl-dev          \
    libcurl4-gnutls-dev \
    libopenmpi-dev      \
    openmpi-bin         \
    libxml2-dev         \
    libfontconfig1-dev  \
    libtool             \
    file                \
    libharfbuzz-dev     \
    libfribidi-dev      \
    libfreetype6-dev    \
    libpng-dev          \
    libtiff5-dev        \
    libjpeg-dev         \
    libuv1-dev          \
    libcairo2-dev       \
    && rm -rf /var/lib/apt/lists/*

RUN Rscript -e "install.packages(c('remotes', 'getopt', 'foreach', 'BiocManager'), repos=Sys.getenv('CRANREPO'))"               \
    && Rscript -e "BiocManager::install('ggtree', update=FALSE, ask=FALSE)"                                                     \
    && Rscript -e "remotes::install_github('emvolz-phylodynamics/treestructure@v0.6.0', dependencies=TRUE, upgrade='never')"    \
    && chmod +x /usr/local/lib/R/site-library/treestructure/tscl


###############################################################################
# Runtime stage
###############################################################################
FROM r-base:4.6.0

LABEL maintainer="anna.zhukova@pasteur.fr"

ENV R_LIBS_SITE=/usr/local/lib/R/site-library
ENV PATH="/usr/local/lib/R/site-library/treestructure:${PATH}"
ENV LC_ALL=C

RUN apt-get update && apt-get install -y --no-install-recommends \
    procps          \
    ca-certificates \
    pandoc          \
    libssl3         \
    libcurl4        \
    libopenmpi40    \
    libxml2         \
    libfontconfig1  \
    libharfbuzz0b   \
    libfribidi0     \
    libfreetype6    \
    libpng16-16     \
    libtiff6        \
    libjpeg62-turbo \
    libuv1          \
    libcairo2       \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir -p /pasteur /work

COPY --from=builder /usr/local/lib/R/site-library /usr/local/lib/R/site-library

WORKDIR /work

CMD ["/bin/bash"]

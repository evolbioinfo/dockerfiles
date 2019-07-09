# R with packages:
# - getopt
# - DBI
# - SRAdb (Bioconductor)
# https://cran.r-project.org/
# v3.6.0

# base image: r-base:3.6.0
FROM r-base:3.6.0

# File Author / Maintainer
MAINTAINER Anna Zhukova <anna.zhukova@pasteur.fr>

ENV CRANREPO="'http://cran.univ-paris1.fr/'"
ENV RLIBPATH="'/usr/local/lib/R/site-library/'"

RUN apt-get update && apt-get install -y procps libssl-dev libcurl4-gnutls-dev curl libopenmpi-dev openmpi-bin openmpi-doc libxml2-dev \
  && R -e "install.packages('DBI',lib=$RLIBPATH,repo=$CRANREPO)" \
  && R -e "install.packages('getopt',lib=$RLIBPATH,repo=$CRANREPO)" \
  && R -e "install.packages('BiocManager',lib=$RLIBPATH,repo=$CRANREPO); BiocManager::install(c('SRAdb'))" \
  && apt-get autoremove -y                                            \
  && apt-get clean                                                    \
  && rm -rf /var/lib/apt/lists/*                                      \
  && mkdir /pasteur

ENTRYPOINT ["/usr/bin/R"]

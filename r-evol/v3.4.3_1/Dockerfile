# R with packages:
# - ape
# - big.phylo
# - ggplot2
# https://cran.r-project.org/
# v3.4.3

# base image: r-base:3.4.3
FROM r-base:3.4.3

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

ENV CRANREPO="'http://cran.univ-paris1.fr/'"
ENV RLIBPATH="'/usr/local/lib/R/site-library/'"

RUN apt-get update && apt-get install -y procps libssl-dev libcurl4-gnutls-dev curl git libopenmpi-dev openmpi-bin openmpi-doc \
  && R -e "install.packages('ape',lib=$RLIBPATH,repo=$CRANREPO)"      \
  && R -e "install.packages('ggplot2',lib=$RLIBPATH,repo=$CRANREPO)"  \
  && R -e "install.packages('devtools',lib=$RLIBPATH,repo=$CRANREPO)" \
  && R -e "library(devtools); with_libpaths(new = $RLIBPATH, install_github(\"olli0601/big.phylo\"))" \
  && apt-get remove curl git -y                                       \
  && apt-get autoremove -y                                            \
  && apt-get clean                                                    \
  && rm -rf /var/lib/apt/lists/*                                      \
  && mkdir /pasteur

ENTRYPOINT ["/usr/bin/R"]

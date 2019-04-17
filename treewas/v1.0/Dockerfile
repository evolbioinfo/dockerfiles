# treeWAS
# https://github.com/caitiecollins/treeWAS
# v1.0

# base image: r-base:3.4.3
FROM r-base:3.4.3

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

ENV CRANREPO="'http://cran.univ-paris1.fr/'"
ENV RLIBPATH="'/usr/local/lib/R/site-library/'"

RUN apt-get update && apt-get install -y procps libmagick++-dev libudunits2-dev libgdal-dev libxml2-dev libssl-dev libcurl4-gnutls-dev \
    && R -e "install.packages('glue',lib=$RLIBPATH,repo=$CRANREPO)" \
    && R -e "install.packages('devtools',lib=$RLIBPATH,repo=$CRANREPO, dep=TRUE)" \
    && R -e "library(devtools); install_github('caitiecollins/treeWAS', build_vignettes = TRUE)" \
    && mkdir /pasteur

ENTRYPOINT ["/usr/bin/R"]

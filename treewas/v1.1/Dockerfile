# treeWAS
# https://github.com/caitiecollins/treeWAS
# v1.1

# base image: r-base:4.5.3
FROM r-base:4.5.3

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

ENV CRANREPO="'https://mirror.ibcp.fr/pub/CRAN/'"
ENV RLIBPATH="'/usr/local/lib/R/site-library/'"

RUN apt-get update                                                                               \
    && apt-get install -y procps libmagick++-dev libbz2-dev liblzma-dev libudunits2-dev libgdal-dev libxml2-dev         \
                          libssl-dev libcurl4-gnutls-dev libgl-dev libglu-dev perl               \
                          libspreadsheet-parseexcel-simple-perl libspreadsheet-parsexlsx-perl    \
                          xlsx2csv catdoc csvkit                                                 \
			  libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
    && apt-get install -y libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev pandoc \
    && R -e "install.packages('devtools',lib=$RLIBPATH,repo=$CRANREPO, dep=TRUE)"                \
    && R -e "library(devtools); install_github('hadley/pryr', build_vignettes = TRUE)"           \
    && R -e "install.packages('glue',lib=$RLIBPATH,repo=$CRANREPO)"                              \
    && R -e "install.packages('gdata',lib=$RLIBPATH,repo=$CRANREPO, dep=TRUE)"                   \
    && R -e "install.packages('phangorn',lib=$RLIBPATH,repo=$CRANREPO, dep=TRUE)"                \
    && R -e "library(devtools); install_github('caitiecollins/treeWAS', build_vignettes = TRUE)" \
    && mkdir /pasteur

ENTRYPOINT ["/usr/bin/R"]

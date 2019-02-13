# script table2itol.R
# Generating iToL annotations from Spreadsheet or CSV files 
# https://github.com/mgoeker/table2itol/
# No version tag - commit e7a38d7545867e0a9ff71781ab39f07a8840c319
# https://github.com/mgoeker/table2itol/tree/e7a38d7545867e0a9ff71781ab39f07a8840c319
# table2itol.R requires R with packages:
# - optparse 
# - plotrix
# - readxl
# - readODS
# - yaml
# https://cran.r-project.org/
# v3.4.3

# base image: r-base:3.4.3
FROM r-base:3.4.3

# File Author / Maintainer
MAINTAINER Mathieu Moslonka-Lefebvre <mathieu.moslonka-lefebvre@pasteur.fr>

ENV CRANREPO="'http://cran.univ-paris1.fr/'"
ENV RLIBPATH="'/usr/local/lib/R/site-library/'"

RUN apt-get update && apt-get install -y procps libxml2 libxml2-dev   \
  && R -e "install.packages('optparse',lib=$RLIBPATH,repo=$CRANREPO)"  \
  && R -e "install.packages('plotrix',lib=$RLIBPATH,repo=$CRANREPO)"     \
  && R -e "install.packages('readxl',lib=$RLIBPATH,repo=$CRANREPO)" \
  && R -e "install.packages('readODS',lib=$RLIBPATH,repo=$CRANREPO)"    \
  && R -e "install.packages('yaml',lib=$RLIBPATH,repo=$CRANREPO)"    \
  && cd /usr/local/ \
  && wget https://github.com/mgoeker/table2itol/archive/e7a38d7545867e0a9ff71781ab39f07a8840c319.zip \
  && unzip e7a38d7545867e0a9ff71781ab39f07a8840c319.zip \
  && cd table2itol-e7a38d7545867e0a9ff71781ab39f07a8840c319 \
  && mv table2itol.R /usr/local/bin/table2itol.R \
  && ln -s /bin/Rscript /usr/local/bin/Rscript \
  && cd /usr/local/ \
  && rm -rf table2itol-e7a38d7545867e0a9ff71781ab39f07a8840c319 \
  && rm -rf e7a38d7545867e0a9ff71781ab39f07a8840c319.zip \
  && apt-get autoremove -y                                            \
  && apt-get clean                                                    \
  && rm -rf /var/lib/apt/lists/*                                      \
  && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/table2itol.R"]

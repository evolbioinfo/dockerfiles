# treedater
# 
# commit 89a0df0  from Jan 16
#
# https://github.com/emvolz/treedater
#
# E.M. Volz and Frost, S.D.W. (2017) Scalable relaxed clock phylogenetic dating. Virus Evolution.

# base image: ubuntu
FROM ubuntu:xenial

# File Author / Maintainer
MAINTAINER Anna Zhukova <anna.zhukova@pasteur.fr>

RUN apt-get update --fix-missing \
    && apt-get install -y wget gcc make bc unzip r-base libssl-dev libcurl4-openssl-dev libssl-dev\
    && apt-get autoremove

ENV CRANREPO="'http://cran.univ-paris1.fr/'"
ENV RLIBPATH="'/usr/local/lib/R/site-library/'"

RUN R -e "install.packages('devtools',lib=$RLIBPATH,repo=$CRANREPO)"  \
   && R -e "install.packages('getopt',lib=$RLIBPATH,repo=$CRANREPO)"  \
  && R -e "library(devtools); install_github( 'emvolz/treedater')"  \
  && apt-get autoremove -y                                            \
  && apt-get clean                                                    \
  && rm -rf /var/lib/apt/lists/*                                      \
  && mkdir /pasteur

ENTRYPOINT ["/usr/local/lib/R/site-library/treedater/tdcl"]

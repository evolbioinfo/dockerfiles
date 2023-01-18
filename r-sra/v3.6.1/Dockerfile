# R with packages:
# - getopt
# - DBI
# - SRAdb (Bioconductor)
# https://cran.r-project.org/
# v3.6.1

# base image: ubuntu
FROM ubuntu:xenial

# File Author / Maintainer
MAINTAINER Anna Zhukova <anna.zhukova@pasteur.fr>

ENV CRANREPO="'http://cran.univ-paris1.fr/'"
ENV RLIBPATH="'/usr/local/lib/R/site-library/'"

RUN apt-get autoremove -y \
    && apt-get update -y \
    && apt-get upgrade -y \
	&& apt-get install -y \
		wget \
		gcc \
		libcurl4-openssl-dev \
		curl \
		libxml2-dev \
		libssl-dev \
		apt-transport-https

## Use Debian unstable via pinning -- new style via APT::Default-Release
RUN echo "deb https://cloud.r-project.org/bin/linux/ubuntu xenial-cran35/" > /etc/apt/sources.list.d/debian-unstable.list

## Now install R
RUN apt-get update -y \
	&& apt-get install -y --allow-unauthenticated r-base-dev

RUN R -e "install.packages('DBI',lib=$RLIBPATH,repo=$CRANREPO)" \
    && R -e "install.packages('getopt',lib=$RLIBPATH,repo=$CRANREPO)" \
    && R -e "install.packages('BiocManager',lib=$RLIBPATH,repo=$CRANREPO); BiocManager::install(c('SRAdb'),lib=$RLIBPATH,repo=$CRANREPO)" \
    && apt-get autoremove -y                                            \
    && apt-get clean                                                    \
    && rm -rf /var/lib/apt/lists/*                                      \
    && mkdir /pasteur

ENTRYPOINT ["/usr/bin/R"]

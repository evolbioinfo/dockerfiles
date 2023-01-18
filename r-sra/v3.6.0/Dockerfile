# R with packages:
# - getopt
# - DBI
# - SRAdb (Bioconductor)
# https://cran.r-project.org/
# v3.6.0

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
		ed \
		less \
		locales \
		vim-tiny \
		wget \
		ca-certificates \
		fonts-texgyre \
		gcc \
		libcurl4-gnutls-dev \
		curl \
		libxml2-dev \
		apt-transport-https \
		ca-certificates \
	&& rm -rf /var/lib/apt/lists/*

## Configure default locale, see https://github.com/rocker-org/rocker/issues/19
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
	&& locale-gen en_US.utf8 \
	&& /usr/sbin/update-locale LANG=en_US.UTF-8

ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8

## Use Debian unstable via pinning -- new style via APT::Default-Release
RUN echo "deb https://cloud.r-project.org/bin/linux/ubuntu xenial-cran35/" > /etc/apt/sources.list.d/debian-unstable.list

RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9

## Now install R
RUN apt-get update -y --allow-unauthenticated \
	&& apt-get dist-upgrade -y --allow-unauthenticated \
	&& apt-get install -y r-base

RUN R -e "install.packages('DBI',lib=$RLIBPATH,repo=$CRANREPO)" \
    && R -e "install.packages('getopt',lib=$RLIBPATH,repo=$CRANREPO)" \
    && R -e "install.packages('BiocManager',lib=$RLIBPATH,repo=$CRANREPO); BiocManager::install(c('SRAdb'))" \
    && apt-get autoremove -y                                            \
    && apt-get clean                                                    \
    && rm -rf /var/lib/apt/lists/*                                      \
    && mkdir /pasteur

ENTRYPOINT ["/usr/bin/R"]

# treestructure
#
# commit a831a66 from Dec 11 2019
#
# https://github.com/emvolz/treestructure
#
# Erik Volz, Carsten Wiuf, Yonatan Grad, Simon Frost, Ann Dennis, Xavier Didelot, bioRxiv 704528;
# doi: https://doi.org/10.1101/704528

# base image: ubuntu
FROM ubuntu:xenial

# File Author / Maintainer
MAINTAINER Anna Zhukova <anna.zhukova@pasteur.fr>

RUN apt-get update --fix-missing \
    && apt-get install -y wget gcc make bc unzip r-base libssl-dev libcurl4-openssl-dev libssl-dev libxml2-dev\
    && apt-get autoremove

ENV CRANREPO="'http://cran.univ-paris1.fr/'"
ENV RLIBPATH="'/usr/local/lib/R/site-library/'"

RUN R -e "install.packages('devtools',lib=$RLIBPATH,repo=$CRANREPO)"  \
   && R -e "install.packages('getopt',lib=$RLIBPATH,repo=$CRANREPO)"  \
   && R -e "install.packages('foreach',lib=$RLIBPATH,repo=$CRANREPO)"  \
  && R -e "library(devtools); install_github( 'emvolz-phylodynamics/treestructure')"  \
  && apt-get autoremove -y                                            \
  && apt-get clean                                                    \
  && rm -rf /var/lib/apt/lists/*                                      \
  && mkdir /pasteur


# create a non-root user named evolbioinfo,
# give them the password "evolbioinfo" put them in the sudo group
RUN useradd -d /home/evolbioinfo -m -s /bin/bash evolbioinfo && echo "evolbioinfo:evolbioinfo" | chpasswd && adduser evolbioinfo sudo

# start working in the "evolbioinfo" home directory
WORKDIR /home/evolbioinfo

# Make the files owned by evolbioinfo
RUN chown -R evolbioinfo:evolbioinfo /home/evolbioinfo

# Switch to your new user in the docker image
USER evolbioinfo

RUN echo "alias tdcl=/usr/local/lib/R/site-library/treestructure/tscl" >> ~/.bashrc

ENTRYPOINT ["/bin/bash"]

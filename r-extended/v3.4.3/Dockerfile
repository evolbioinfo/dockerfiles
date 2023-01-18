# R with packages:
# - ggplot2
# - plyr
# - reshape2
# - dplyr, tidyr
# https://cran.r-project.org/
# v3.4.3

# base image: r-base:3.4.3
FROM r-base:3.4.3

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

ENV CRANREPO="'http://cran.univ-paris1.fr/'"
ENV RLIBPATH="'/usr/local/lib/R/site-library/'"

RUN apt-get update && apt-get install -y procps                       \
  && R -e "install.packages('ggplot2',lib=$RLIBPATH,repo=$CRANREPO)"  \
  && R -e "install.packages('plyr',lib=$RLIBPATH,repo=$CRANREPO)"     \
  && R -e "install.packages('reshape2',lib=$RLIBPATH,repo=$CRANREPO)" \
  && R -e "install.packages('dplyr',lib=$RLIBPATH,repo=$CRANREPO)"    \
  && R -e "install.packages('tidyr',lib=$RLIBPATH,repo=$CRANREPO)"    \
  && apt-get autoremove -y                                            \
  && apt-get clean                                                    \
  && rm -rf /var/lib/apt/lists/*                                      \
  && mkdir /pasteur

ENTRYPOINT ["/usr/bin/R"]

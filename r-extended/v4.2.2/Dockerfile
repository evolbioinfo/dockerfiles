# R with packages:
# - ggplot2
# - plyr
# - reshape2
# - dplyr, tidyr
# - lubridate
# - data.table
# - rmarkdown
# - plotly
# https://cran.r-project.org/
# v4.2.2

# base image: r-base:4.2.2
FROM r-base:4.2.2

# File Author / Maintainer
LABEL maintainer=frederic.lemoine@pasteur.fr

ENV CRANREPO="'https://pbil.univ-lyon1.fr/CRAN/'"
ENV RLIBPATH="'/usr/local/lib/R/site-library/'"

RUN apt-get update && apt-get install -y procps libssl-dev libcurl4-openssl-dev \
  && R -e "install.packages('ggplot2',lib=$RLIBPATH,repo=$CRANREPO)"  \
  && R -e "install.packages('plyr',lib=$RLIBPATH,repo=$CRANREPO)"     \
  && R -e "install.packages('reshape2',lib=$RLIBPATH,repo=$CRANREPO)" \
  && R -e "install.packages('dplyr',lib=$RLIBPATH,repo=$CRANREPO)"    \
  && R -e "install.packages('tidyr',lib=$RLIBPATH,repo=$CRANREPO)"    \
  && R -e "install.packages('lubridate',lib=$RLIBPATH,repo=$CRANREPO)" \
  && R -e "install.packages('data.table',lib=$RLIBPATH,repo=$CRANREPO)" \
  && R -e "install.packages('plotly',lib=$RLIBPATH,repo=$CRANREPO)"   \
  && R -e "install.packages('rmarkdown',lib=$RLIBPATH,repo=$CRANREPO)"\
  && apt-get autoremove -y                                            \
  && apt-get clean                                                    \
  && rm -rf /var/lib/apt/lists/*                                      \
  && mkdir /pasteur

ENTRYPOINT ["/usr/bin/R"]

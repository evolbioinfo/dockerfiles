# R with packages:
# - ape
# - big.phylo
# - ggplot2
# - plyr
# - dplyr
# - tidyr
# - reshape2
# - lubridate
# - Phylogram
# - Admixtools
# https://cran.r-project.org/
# v4.1.0

# base image: r-base:4.1.0
FROM r-base:4.1.0

# File Author / Maintainer
LABEL maintainer=frederic.lemoine@pasteur.fr

ENV CRANREPO="'http://cran.univ-paris1.fr/'"
ENV RLIBPATH="'/usr/local/lib/R/site-library/'"

RUN apt-get update && apt-get install -y procps libssl-dev            \
              libcurl4-gnutls-dev curl git libopenmpi-dev             \
	      openmpi-bin openmpi-doc libxml2-dev apt-transport-https \
  && R -e "install.packages('devtools',lib=$RLIBPATH,repo=$CRANREPO)" \
  && R -e "library(devtools); withr::with_libpaths(new = $RLIBPATH, install_github(\"uqrmaie1/admixtools\",INSTALL_opts=\"--no-staged-install\"))" \
  && R -e "install.packages('ape',lib=$RLIBPATH,repo=$CRANREPO)"      \
  && R -e "install.packages('treedater',lib=$RLIBPATH,repo=$CRANREPO)"\
  && R -e "install.packages('phylogram',lib=$RLIBPATH,repo=$CRANREPO)"\
  && R -e "install.packages('ggplot2',lib=$RLIBPATH,repo=$CRANREPO)"  \
  && R -e "install.packages('plyr',lib=$RLIBPATH,repo=$CRANREPO)"     \
  && R -e "install.packages('reshape2',lib=$RLIBPATH,repo=$CRANREPO)" \
  && R -e "install.packages('dplyr',lib=$RLIBPATH,repo=$CRANREPO)"    \
  && R -e "install.packages('tidyr',lib=$RLIBPATH,repo=$CRANREPO)"    \
  && R -e "install.packages('lubridate',lib=$RLIBPATH,repo=$CRANREPO)"\
  && R -e "library(devtools); withr::with_libpaths(new = $RLIBPATH, install_github(\"olli0601/big.phylo\",INSTALL_opts=\"--no-staged-install\"))" \
  && apt-get autoremove -y                                            \
  && apt-get clean                                                    \
  && rm -rf /var/lib/apt/lists/*                                      \
  && mkdir /pasteur

ENTRYPOINT ["/usr/bin/R"]

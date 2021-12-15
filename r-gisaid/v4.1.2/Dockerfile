# R with packages:
# lubridate
# devtools
# GISAIDR
# v4.1.2

# base image: r-base:4.1.2
FROM r-base:4.1.2

# File Author / Maintainer
LABEL maintainer=julien.fumey@pasteur.fr

ENV CRANREPO="'http://cran.univ-paris1.fr/'"
ENV RLIBPATH="'/usr/local/lib/R/site-library'"

RUN apt-get update && apt-get install -y procps \
    && R -e "install.packages('lubridate', lib=$RLIBPATH, repo=$CRANREPO)"  \
    && R -e "install.packages('devtools', lib=$RLIBPATH, repo=$CRANREPO)"   \
    && R -e "library(devtools); withr::with_libpaths(new = $RLIBPATH, install_github(\"julienfumey/GISAIDR\",INSTALL_opts=\"--no-staged-install\"))" \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/bin/R"]
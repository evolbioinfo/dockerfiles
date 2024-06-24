# R with packages:
# - ggplot2
# - plyr
# - reshape2
# - dplyr, tidyr
# - lubridate
# - data.table
# - rmarkdown
# - plotly
# - ggthemes
# https://cran.r-project.org/
# v4.3.3

# base image: r-base:4.3.3
FROM r-base:4.3.3

# File Author / Maintainer
LABEL maintainer=frederic.lemoine@pasteur.fr

ENV CRANREPO="https://pbil.univ-lyon1.fr/CRAN/"
ENV RLIBPATH="/usr/local/lib/R/site-library/"

COPY install.R /
COPY packages.txt /

RUN apt-get update && apt-get install -y procps libssl-dev  pandoc \
    pandoc pandoc-data pandoc-citeproc-preamble pandoc-plantuml-filter libghc-pandoc-sidenote-dev \
    libghc-pandoc-sidenote-doc libghc-pandoc-sidenote-prof pandoc-sidenote \
    ghostscript lmodern qpdf texinfo texlive-fonts-extra texlive-fonts-recommended \
    texlive-latex-extra texlive-latex-recommended texlive-luatex texlive-plain-generic libfreetype6-dev \
    libfontconfig1-dev libxml2-dev texlive-science texlive-xetex libharfbuzz-dev libfribidi-dev \
    libcurl4-gnutls-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev libmagick++-dev libudunits2-dev cargo libgdal-dev libnlopt-dev libtool m4 automake

RUN apt-get install -y libv8-dev

RUN export MAKE="make -j15" \
  && Rscript /install.R $CRANREPO $RLIBPATH /packages.txt  \
  && R -e "devtools::install_github('NightingaleHealth/ggforestplot')" \
  && apt-get autoremove -y                                            \
  && apt-get clean                                                    \
  && rm -rf /var/lib/apt/lists/*                                      \
  && mkdir /pasteur

ENTRYPOINT ["/usr/bin/R"]


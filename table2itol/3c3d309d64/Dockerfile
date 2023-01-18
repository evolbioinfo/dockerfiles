# script table2itol.R
# Generating iToL annotations from Spreadsheet or CSV files 
# https://github.com/mgoeker/table2itol/
# No version tag - commit 3c3d309d64c515797b63835f330b41e1294a1f7f
# https://github.com/mgoeker/table2itol/tree/3c3d309d64c515797b63835f330b41e1294a1f7f
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

RUN apt-get update && apt-get install -y procps                       \
  && R -e "install.packages('optparse',lib=$RLIBPATH,repo=$CRANREPO)"  \
  && R -e "install.packages('plotrix',lib=$RLIBPATH,repo=$CRANREPO)"     \
  && R -e "install.packages('readxl',lib=$RLIBPATH,repo=$CRANREPO)" \
  && R -e "install.packages('readODS',lib=$RLIBPATH,repo=$CRANREPO)"    \
  && R -e "install.packages('yaml',lib=$RLIBPATH,repo=$CRANREPO)"    \
  && cd /usr/local/ \
  && wget https://github.com/mgoeker/table2itol/archive/3c3d309d64c515797b63835f330b41e1294a1f7f.zip \
  && unzip 3c3d309d64c515797b63835f330b41e1294a1f7f.zip \
  && cd table2itol-3c3d309d64c515797b63835f330b41e1294a1f7f \
  && mv table2itol.R /usr/local/bin/table2itol.R \
  && ln -s /bin/Rscript /usr/local/bin/Rscript \
  && cd /usr/local/ \
  && rm -rf table2itol-3c3d309d64c515797b63835f330b41e1294a1f7f \
  && rm -rf 3c3d309d64c515797b63835f330b41e1294a1f7f.zip \
  && apt-get autoremove -y                                            \
  && apt-get clean                                                    \
  && rm -rf /var/lib/apt/lists/*                                      \
  && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/table2itol.R"]

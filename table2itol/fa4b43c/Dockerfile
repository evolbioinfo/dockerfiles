# script table2itol.R
# Generating iToL annotations from Spreadsheet or CSV files 
# https://github.com/mgoeker/table2itol/
# No version tag - commit fa4b43c16c8bc437598c0d299fba96bc78596917
# https://github.com/mgoeker/table2itol/tree/fa4b43c16c8bc437598c0d299fba96bc78596917
# table2itol.R requires R with packages:
# - optparse 
# - plotrix
# - readxl
# - readODS
# - yaml
# https://cran.r-project.org/
# v4.0.2

# base image: r-base:4.0.2

FROM r-base:4.0.2

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

ENV CRANREPO="'https://pbil.univ-lyon1.fr/CRAN/'"
ENV RLIBPATH="'/usr/local/lib/R/site-library/'"

RUN apt-get update && apt-get install -y apt-utils && apt-get install -y procps libxml2 libxml2-dev \
  && R -e "install.packages('optparse',lib=$RLIBPATH,repo=$CRANREPO)"  \
  && R -e "install.packages('plotrix',lib=$RLIBPATH,repo=$CRANREPO)"     \
  && R -e "install.packages('readxl',lib=$RLIBPATH,repo=$CRANREPO)" \
  && R -e "install.packages('readODS',lib=$RLIBPATH,repo=$CRANREPO)"    \
  && R -e "install.packages('yaml',lib=$RLIBPATH,repo=$CRANREPO)"    \
  && cd /usr/local/ \
  && wget https://github.com/mgoeker/table2itol/archive/fa4b43c16c8bc437598c0d299fba96bc78596917.zip \
  && unzip fa4b43c16c8bc437598c0d299fba96bc78596917.zip \
  && cd table2itol-fa4b43c16c8bc437598c0d299fba96bc78596917 \
  && mv table2itol.R /usr/local/bin/table2itol.R \
  && ln -s /bin/Rscript /usr/local/bin/Rscript \
  && cd /usr/local/ \
  && rm -rf table2itol-fa4b43c16c8bc437598c0d299fba96bc78596917 \
  && rm -rf fa4b43c16c8bc437598c0d299fba96bc78596917.zip \
  && apt-get autoremove -y                                            \
  && apt-get clean                                                    \
  && rm -rf /var/lib/apt/lists/*                                      \
  && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/table2itol.R"]

# mapDamage
# https://ginolhac.github.io/mapDamage/
# v2.2.1

# base image: r-base:3.4.3
FROM r-base:4.0.3

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

ENV CRANREPO="'http://cran.univ-paris1.fr/'"
ENV RLIBPATH="'/usr/local/lib/R/site-library/'"

RUN apt-get update \
    && apt-get install -y procps git python3 libgeos-dev \
        file libgsl-dev make gcc libz-dev libudunits2-dev \
        python3-dev python3-pip
RUN pip install pysam
RUN cd /usr/local/ && git clone https://github.com/ginolhac/mapDamage.git && cd mapDamage && git checkout 2.2.1
RUN cd /usr/local/mapDamage && python3 setup.py install
RUN R -e "install.packages('inline',lib=$RLIBPATH,repo=$CRANREPO)"                   
RUN R -e "install.packages('gam',lib=$RLIBPATH,repo=$CRANREPO)"                             
RUN R -e "install.packages('Rcpp',lib=$RLIBPATH,repo=$CRANREPO, dep=TRUE)"                   
RUN R -e "install.packages('ggplot2',lib=$RLIBPATH,repo=$CRANREPO, dep=TRUE)"
RUN R -e "install.packages('RcppGSL',lib=$RLIBPATH,repo=$CRANREPO, dep=TRUE)"
RUN mkdir /pasteur

ENTRYPOINT ["/usr/bin/R"]

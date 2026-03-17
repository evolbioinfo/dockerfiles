# mapDamage
# https://ginolhac.github.io/mapDamage/
# https://github.com/ginolhac/mapDamage/
# v2.2.3

# base image: r-base:3.4.3
FROM r-base:4.5.3

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

ENV CRANREPO="'http://cran.univ-paris1.fr/'"
ENV RLIBPATH="'/usr/local/lib/R/site-library/'"

RUN apt-get update \
    && apt-get  install -y procps git python3 libgeos-dev \
        file libgsl-dev make gcc libz-dev libudunits2-dev \
        python3-dev python3-pip libcurl4-gnutls-dev \
	libfontconfig1-dev libxml2-dev libssl-dev libcairo2-dev git \
	libopenmpi-dev openmpi-bin openmpi-doc wget libbz2-dev zlib1g \
	zlib1g-dev liblzma5 liblzma-dev ncurses-base libncurses-dev bzip2 \
	libgit2-dev seqtk python3-pysam python3-setuptools \
    && cd /usr/local/ \
    && R -e "install.packages('inline',lib=$RLIBPATH,repo=$CRANREPO)" \ 
    && R -e "install.packages('gam',lib=$RLIBPATH,repo=$CRANREPO)" \
    && R -e "install.packages('Rcpp',lib=$RLIBPATH,repo=$CRANREPO, dep=TRUE)" \
    && R -e "install.packages('ggplot2',lib=$RLIBPATH,repo=$CRANREPO, dep=TRUE)" \
    && R -e "install.packages('RcppGSL',lib=$RLIBPATH,repo=$CRANREPO, dep=TRUE)" \
    && cd /usr/local/ && git clone https://github.com/ginolhac/mapDamage.git && cd mapDamage && git checkout 2.2.3 && python3 setup.py install \
    && mkdir /pasteur

ENTRYPOINT ["/usr/bin/R"]

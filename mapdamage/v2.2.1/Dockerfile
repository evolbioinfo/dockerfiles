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
        python3-dev python3-pip libcurl4-openssl-dev \
	libfontconfig1-dev libxml2-dev libssl-dev libcairo2-dev curl git \
	libopenmpi-dev openmpi-bin openmpi-doc wget libbz2-dev zlib1g \
	zlib1g-dev liblzma5 liblzma-dev libncurses5 libncurses5-dev bzip2 libgit2-dev \
    && cd /usr/local/ \
    && wget https://github.com/lh3/seqtk/archive/v1.3.tar.gz \
    && tar -xzvf v1.3.tar.gz \
    && rm v1.3.tar.gz \
    && cd seqtk-1.3 / \
    && make \
    && make install \
    && cd /usr/local \
    && rm -rf /usr/local/seqtk-1.3 \
    && pip install pysam \
    && R -e "install.packages('inline',lib=$RLIBPATH,repo=$CRANREPO)" \ 
    && R -e "install.packages('gam',lib=$RLIBPATH,repo=$CRANREPO)" \
    && R -e "install.packages('Rcpp',lib=$RLIBPATH,repo=$CRANREPO, dep=TRUE)" \
    && R -e "install.packages('ggplot2',lib=$RLIBPATH,repo=$CRANREPO, dep=TRUE)" \
    && R -e "install.packages('RcppGSL',lib=$RLIBPATH,repo=$CRANREPO, dep=TRUE)" \
    && cd /usr/local/ && git clone https://github.com/ginolhac/mapDamage.git && cd mapDamage && git checkout 2.2.1 && python3 setup.py install \
    && mkdir /pasteur

ENTRYPOINT ["/usr/bin/R"]

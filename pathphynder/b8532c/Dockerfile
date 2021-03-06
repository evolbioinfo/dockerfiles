# pathPhynder:
# https://github.com/ruidlpm/pathPhynder
# R-base : v4.0.2

# base image: r-base:4.0.2
FROM r-base:4.0.2

# File Author / Maintainer
LABEL maintainer=alemeur@u-psud.fr

ENV CRANREPO="'http://cran.univ-paris1.fr/'"
ENV RLIBPATH="'/usr/local/lib/R/site-library/'"
ENV VERSION="b8532c07b7292812ae4ddc6246c68eb8e2534b7f"
ENV PHYNDERVERSION="48b184f76c67985ef019ce199e1a1a7dafdc51a5"
ENV HTSLIBVERSION="1.12"
ENV SAMTOOLSVERSION="1.12"

COPY pathPhynder.sh /usr/local/bin/pathPhynder
RUN apt-get update && apt-get install -y git libbz2-dev zlib1g zlib1g-dev liblzma5 liblzma-dev \
    	    	      	      	         libncurses5 libncurses5-dev python3 libcurl4-gnutls-dev procps \
					 gcc wget make libbz2-dev bzip2 \
	&& R -e "install.packages('maps',lib=$RLIBPATH,repo=$CRANREPO)" \
	&& R -e "install.packages('phytools',lib=$RLIBPATH,repo=$CRANREPO)" \
	&& R -e "install.packages('optparse',lib=$RLIBPATH,repo=$CRANREPO)" \
	&& R -e "install.packages('scales',lib=$RLIBPATH,repo=$CRANREPO)" \
	&& cd /usr/local/ \
	&& git clone https://github.com/ruidlpm/pathPhynder.git \
	&& cd pathPhynder \
	&& git checkout $VERSION \
	&& echo 'alias pathPhynder="Rscript /usr/local/pathPhynder/pathPhynder.R"' >> /etc/profile \
	&& cd /usr/local/ \
	&& git clone https://github.com/samtools/htslib.git \
	&& cd htslib \
	&& git checkout $HTSLIBVERSION \
	&& git submodule update --init --recursive \
	&& make lib-static htslib_static.mk \
        && cd /usr/local/ \
        && wget -O samtools.tar.bz2 https://github.com/samtools/samtools/releases/download/${SAMTOOLSVERSION}/samtools-${SAMTOOLSVERSION}.tar.bz2 \
        && tar -xjvf samtools.tar.bz2 \
        && rm -rf samtools.tar.bz2 \
        && cd samtools-${SAMTOOLSVERSION} \
        && ./configure \
        && make \
        && make install \
        && cd /usr/local \
        && rm -rf /usr/local/samtools-${SAMTOOLSVERSION} \
	&& cd /usr/local/ \
        && git clone https://github.com/richarddurbin/phynder.git \
	&& cd phynder \
	&& git checkout $PHYNDERVERSION \
	&& make \
	&& chmod 777 phynder \
	&& mv phynder /usr/local/bin/ \
	&& cd .. \
	&& rm -rf phynder \
	&& chmod +x /usr/local/bin/pathPhynder \
	&& apt-get remove -y git \
	&& apt-get autoremove -y \
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/* \
	&& mkdir /pasteur

ENTRYPOINT ["pathPhynder"]

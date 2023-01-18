# Minion_QC
# https://github.com/roblanf/minion_qc
# v1.4.1

# base image: r-base:3.4.3
FROM r-base:3.4.3

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

ENV CRANREPO="'http://cran.univ-paris1.fr/'"
ENV RLIBPATH="'/usr/local/lib/R/site-library/'"
ENV PATH="/usr/local/minion_qc/:${PATH}"

RUN apt-get update \
    && apt-get install -y git \
    && git clone https://github.com/roblanf/minion_qc.git \
    && cd minion_qc \
    && git checkout v1.4.1 \
    && rm -rf .git \
    && cd .. \
    && cp -r minion_qc /usr/local/ \
    && chmod +x /usr/local/minion_qc/MinIONQC.R /usr/local/minion_qc/tests.R \
    && R -e 'install.packages(c("data.table","futile.logger","ggplot2","optparse","plyr","readr","reshape2","scales","viridis","yaml"))' \
    && apt-get remove -y git \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/minion_qc/MinIONQC.R"]

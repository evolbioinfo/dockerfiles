# Strainline
# https://github.com/HaploKit/Strainline
# Commit 8af032906e75469e77e4b7d7c325a7493c6f5bb5

FROM continuumio/miniconda3:4.10.3

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

ENV VERSION="8af032906e75469e77e4b7d7c325a7493c6f5bb5"

COPY file.patch /

RUN apt-get update --fix-missing \
    && apt-get install -y  git patch \
    && cd / \
    && conda config --add channels defaults \
    && conda config --add channels bioconda \
    && conda config --add channels conda-forge \
    && conda install -c bioconda minimap2==2.24 spoa==4.0.7 samtools==1.14 dazz_db==1.0 daligner==1.0.20200322 metabat2==2.12.1 \
    && wget https://github.com/gt1/daccord/releases/download/0.0.10-release-20170526170720/daccord-0.0.10-release-20170526170720-x86_64-etch-linux-gnu.tar.gz \
    && tar -zvxf daccord-0.0.10-release-20170526170720-x86_64-etch-linux-gnu.tar.gz \
    && rm -rf daccord-0.0.10-release-20170526170720-x86_64-etch-linux-gnu.tar.gz \
    && ln -fs $PWD/daccord-0.0.10-release-20170526170720-x86_64-etch-linux-gnu/bin/daccord /opt/conda/bin/daccord \
    && git clone https://github.com/HaploKit/Strainline.git \
    && cd Strainline \
    && git checkout $VERSION \
    && cd src \
    && patch -f < /file.patch \
    && mv * /opt/conda/bin/ \
    && cd ../../ \
    && rm -rf Strainline

WORKDIR /tmp/

ENTRYPOINT ["/opt/conda/bin/strainline.sh"]

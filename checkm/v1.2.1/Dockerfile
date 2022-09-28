# python3.9 with packages:
# - numpy=1.23.3
# - pysam=0.19.1
# - prodigal=2.6.3
# - hmmer=3.3.2
# - pplacer=v1.1.alpha19
# - matplotlib=3.6.0
# - pysam=0.19.1
# - checkm-genome=1.2.1

FROM ubuntu:22.04

ENV PATH=/miniconda/bin:${PATH}

# Install CheckM
RUN apt-get update --fix-missing \
    && apt-get install -y wget gcc \
    && wget http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && bash Miniconda3-latest-Linux-x86_64.sh -p /miniconda -b \
    && rm Miniconda3-latest-Linux-x86_64.sh \
    && conda install -c conda-forge -y  mamba \
    && mamba install -c anaconda -y python=3.9 \
    && mamba install -y pip \
    && mamba install -c bioconda -y prodigal=2.6.3 hmmer=3.3.2 pplacer=1.1.alpha19 \
    && pip install numpy==1.23.3 matplotlib==3.6.0 pysam==0.19.1 checkm-genome==1.2.1 \
    && mkdir /checkm_data \
    && cd /checkm_data \
    && wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz \
    && tar -xzvf checkm_data_2015_01_16.tar.gz \
    && checkm data setRoot /checkm_data

# The entrypoint runs CheckM
ENTRYPOINT ["checkm"]

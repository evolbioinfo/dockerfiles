# python3.8.2 with packages:
# - numpy=1.20.1
# - pysam=0.16.0.1
# - prodigal=2.6.3
# - hmmer=3.3.2
# - pplacer=v1.1.alpha19
# - matplotlib=3.3.4
# - pysam=0.16.0.1
# - checkm-genome=1.1.3

FROM ubuntu:20.04

ENV PATH=/miniconda/bin:${PATH}

# Install CheckM
RUN apt-get update --fix-missing \
    && apt-get install -y wget gcc \
    && wget http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && bash Miniconda3-latest-Linux-x86_64.sh -p /miniconda -b \
    && rm Miniconda3-latest-Linux-x86_64.sh \
    && conda update -y conda \
    && conda install -c anaconda -y python=3.8.2 \
    && conda install pip \
    && conda install -c bioconda prodigal=2.6.3 hmmer=3.3.2 pplacer=1.1.alpha19 \
    && pip install numpy==1.20.1 matplotlib==3.3.4 pysam==0.16.0.1 checkm-genome==1.1.3 \
    && mkdir /checkm_data \
    && cd /checkm_data \
    && wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz \
    && tar -xzvf checkm_data_2015_01_16.tar.gz \
    && checkm data setRoot /checkm_data

# The entrypoint runs CheckM
ENTRYPOINT ["checkm"]

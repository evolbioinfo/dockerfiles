# python3.9 with packages:
# - numpy=1.23.3
# - pysam=0.19.1
# - prodigal=2.6.3
# - hmmer=3.3.2
# - pplacer=v1.1.alpha19
# - matplotlib=3.6.0
# - pysam=0.19.1
# - checkm-genome=1.2.5

FROM continuumio/miniconda3:4.12.0

# Install CheckM
RUN apt-get update --fix-missing \
    && apt-get install -y wget gcc \
    && conda install -c bioconda -y prodigal=2.6.3 hmmer=3.3.2 pplacer=1.1.alpha19 \
    && pip install numpy==1.23.3 matplotlib==3.6.0 pysam==0.19.1 checkm-genome==1.2.5 \
    && mkdir /checkm_data \
    && cd /checkm_data \
    && wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz \
    && tar -xzvf checkm_data_2015_01_16.tar.gz \
    && checkm data setRoot /checkm_data

# The entrypoint runs CheckM
ENTRYPOINT ["checkm"]

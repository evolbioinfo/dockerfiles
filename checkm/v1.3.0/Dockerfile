# CheckM
# https://github.com/Ecogenomics/CheckM
# v1.3.0
# Note: CheckM 1.3.0 replaces pplacer by gPlace and changes some options and
# output files (e.g. quality_report.tsv, gzip-compressed outputs).
# See https://github.com/Ecogenomics/CheckM/releases/tag/v1.3.0
# The bioconda package brings python, prodigal, hmmer and gtdb-gplace.

FROM condaforge/miniforge3:26.1.0-0

# File Author / Maintainer
LABEL maintainer=frederic.lemoine@pasteur.fr

ENV VERSION="1.3.0"
ENV PATH=/opt/conda/envs/checkm/bin:$PATH

# Install CheckM in its own environment, and its reference data
RUN conda create -n checkm -c conda-forge -c bioconda -y checkm-genome=${VERSION} \
    && conda clean -a -y \
    && mkdir /checkm_data \
    && cd /checkm_data \
    && wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz \
    && tar -xzvf checkm_data_2015_01_16.tar.gz \
    && rm -f checkm_data_2015_01_16.tar.gz \
    && checkm data setRoot /checkm_data \
    && mkdir /pasteur

# The entrypoint runs CheckM
ENTRYPOINT ["checkm"]

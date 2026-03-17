# pathPhynder:
# https://github.com/ruidlpm/pathPhynder
# v1.2.4

FROM condaforge/miniforge3:26.1.0-0

# File Author / Maintainer
LABEL maintainer=frederic.lemoine@pasteur.fr

ENV VERSION="1.2.4"

# Install PathPhynder
RUN  conda install -c bioconda -y pathphynder=${VERSION} samtools=1.20

# The entrypoint runs CheckM
ENTRYPOINT ["pathPhynder"]

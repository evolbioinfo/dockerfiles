# pathPhynder:
# https://github.com/ruidlpm/pathPhynder
# v1.3.1

FROM condaforge/miniforge3:26.1.0-0

# File Author / Maintainer
LABEL maintainer=frederic.lemoine@pasteur.fr

ENV VERSION="1.3.1"
ENV PATH=/opt/conda/envs/pathphynder/bin:$PATH

# Install PathPhynder in its own environment (pathphynder 1.3.1 needs
# python 3.14, which the base environment does not have)
RUN conda create -n pathphynder -c conda-forge -c bioconda -y pathphynder=${VERSION} samtools=1.20 \
    && conda clean -a -y \
    && mkdir /pasteur

# The entrypoint runs pathPhynder
ENTRYPOINT ["pathPhynder"]

# bwa-mem2
# https://github.com/bwa-mem2/bwa-mem2
# https://anaconda.org/bioconda/bwa-mem2
# v2.2.1

# base image: Ubuntu
FROM ubuntu:24.04

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

ENV VERSION=2.2.1
ENV SAMVERSION=1.23.1


# Use a minimal Miniconda base image (available for amd64 and arm64)
FROM condaforge/miniforge3:latest

# Install bwa-mem2 from bioconda
RUN conda install -y -c bioconda -c conda-forge bwa-mem2=${VERSION} samtools=${SAMVERSION} \
    && conda clean --all --yes

# Set a default working directory
WORKDIR /data

# Set the default command
ENTRYPOINT ["bwa-mem2"]
CMD ["--help"]

# civet
# https://github.com/artic-network/polecat
# commit b4a36f34dc523d0c5b1d829b26d9ff6e2a6cc21c

FROM ubuntu:20.04

ENV PATH=/miniconda/envs/polecat/bin:${PATH}

# Install civet
RUN apt-get update --fix-missing \
    && apt-get install -y wget git gcc \
    && wget http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && bash Miniconda3-latest-Linux-x86_64.sh -p /miniconda -b \
    && rm Miniconda3-latest-Linux-x86_64.sh \
    && git clone https://github.com/artic-network/polecat/ \
    && cd polecat \
    && git checkout b4a36f34dc523d0c5b1d829b26d9ff6e2a6cc21c \
    && /miniconda/bin/conda env create -f environment.yml \
    && python setup.py install \
    && cd .. \
    && rm -rf polecat \
    && apt-get remove -y wget git gcc \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

# The entrypoint runs civet
ENTRYPOINT ["polecat"]

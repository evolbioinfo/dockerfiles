# civet
# https://github.com/artic-network/civet/
# version v2.1.2

FROM ubuntu:20.04

ENV PATH=/miniconda/envs/civet/bin:${PATH}

# Install civet
RUN apt-get update --fix-missing \
    && apt-get install -y wget git \
    && wget http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && bash Miniconda3-latest-Linux-x86_64.sh -p /miniconda -b \
    && rm Miniconda3-latest-Linux-x86_64.sh \
    && git clone https://github.com/COG-UK/civet.git \
    && cd civet \
    && git checkout v2.1.2 \
    && /miniconda/bin/conda env create -f environment.yml \
    && python setup.py install \
    && cd .. \
    && rm -rf civet \
    && apt-get remove -y wget git \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

# The entrypoint runs civet
ENTRYPOINT ["civet"]

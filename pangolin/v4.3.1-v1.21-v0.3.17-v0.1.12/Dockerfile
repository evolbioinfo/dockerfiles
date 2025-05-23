# Largely inspired from
#  https://github.com/cov-lineages/pangolin/blob/master/Dockerfile

FROM continuumio/miniconda3:23.5.2-0-alpine
LABEL version="4.3.1" \
      description="Docker image for Pangolin"

ENV VERSION "v4.3.1"
ENV DATAVERSION "v1.21"
ENV SCORPIOVERSION "v0.3.17"
ENV CONSTELLATIONVERSION "v0.1.12"

# Install git for pangolin
RUN apk update \
    && apk add git bash libarchive \
    && pip install --upgrade pip \
    && pip uninstall --yes pyOpenSSL \
    && pip install --upgrade cryptography \
    && pip install pyOpenSSL \
    && conda update conda \
    && conda install mamba -c conda-forge

RUN mamba config --set safety_checks disabled

RUN cd / \
    && git clone https://github.com/cov-lineages/pangolin.git \
    && cd pangolin \
    && git checkout $VERSION \
    && cd .. \
    && sed -i 's/snakemake-minimal=/snakemake-minimal>=/' /pangolin/environment.yml \
    && sed -i 's/pip=19.3.1/pip>=19.3.1/' /pangolin/environment.yml \
    && sed -i 's/name: pangolin/name: base/' /pangolin/environment.yml \
    && sed -i s=git+https://github.com/cov-lineages/pangolin-data.git=git+https://github.com/cov-lineages/pangolin-data.git@$DATAVERSION= /pangolin/environment.yml \
    && sed -i s=git+https://github.com/cov-lineages/scorpio.git=git+https://github.com/cov-lineages/scorpio.git@$SCORPIOVERSION= /pangolin/environment.yml \
    && sed -i s=git+https://github.com/cov-lineages/constellations.git=git+https://github.com/cov-lineages/constellations.git@$CONSTELLATIONVERSION= /pangolin/environment.yml
    
#RUN mamba install -y wheel certifi six pycparser jsonpointer charset-normalizer jsonpatch --force-reinstall -c conda-forge \
#    && mamba env update  --quiet -f /pangolin/environment.yml && conda clean -a \
#    && cd /pangolin \
#    && pip install . && rm -rf /root/.cache/pip \
#    && pangolin --version &> /pangolin-version.txt
    
WORKDIR /tmp/

ENTRYPOINT ["pangolin"]

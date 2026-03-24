# Largely inspired from
# https://github.com/cov-lineages/pangolin/blob/master/Dockerfile

FROM condaforge/miniforge3:26.1.0-0

ARG DEBIAN_FRONTEND=noninteractive
LABEL version="4.4" \
      description="Docker image for Pangolin"

ENV VERSION=v4.4
ENV DATAVERSION=v1.38
ENV SCORPIOVERSION=v0.3.19
ENV CONSTELLATIONVERSION=v0.1.12
ENV CONDA_ENV_NAME=pangolin
ENV PATH=/opt/conda/envs/pangolin/bin:$PATH

COPY pangolin.patch /

# Install git for pangolin
RUN apt-get update -y \
    && apt-get install -y git bash libarchive-dev \
    && pip install --upgrade pip \
    && pip uninstall --yes pyOpenSSL \
    && pip install --upgrade cryptography \
    && pip install pyOpenSSL \
    && conda clean -afy

RUN cd / \
    && git clone https://github.com/cov-lineages/pangolin.git \
    && cd pangolin \
    && git checkout $VERSION \
    && git apply /pangolin.patch \
    && cd .. \
    && sed -i s=git+https://github.com/cov-lineages/pangolin-data.git=git+https://github.com/cov-lineages/pangolin-data.git@$DATAVERSION= /pangolin/environment.yml \
    && sed -i s=git+https://github.com/cov-lineages/scorpio.git=git+https://github.com/cov-lineages/scorpio.git@$SCORPIOVERSION= /pangolin/environment.yml \
    && sed -i s=git+https://github.com/cov-lineages/constellations.git=git+https://github.com/cov-lineages/constellations.git@$CONSTELLATIONVERSION= /pangolin/environment.yml \
    && sed -i 's/^  - python>=3.10$/  - python=3.11/' /pangolin/environment.yml \
    && sed -i 's/^  - pip>=19.3.1$/  - pip/' /pangolin/environment.yml \
    && sed -i '/^  - pip:/,$d' /pangolin/environment.yml \
    && conda env create --quiet -n $CONDA_ENV_NAME -f /pangolin/environment.yml \
    && conda clean -afy \
    && cd /pangolin \
    && conda run --no-capture-output -n $CONDA_ENV_NAME pip install --no-cache-dir 'setuptools<81' wheel \
    && conda run --no-capture-output -n $CONDA_ENV_NAME pip install --no-cache-dir --no-build-isolation git+https://github.com/cov-lineages/pangolin-data.git@$DATAVERSION \
    && conda run --no-capture-output -n $CONDA_ENV_NAME pip install --no-cache-dir --no-build-isolation git+https://github.com/cov-lineages/scorpio.git@$SCORPIOVERSION \
    && conda run --no-capture-output -n $CONDA_ENV_NAME pip install --no-cache-dir --no-build-isolation git+https://github.com/cov-lineages/constellations.git@$CONSTELLATIONVERSION \
    && conda run --no-capture-output -n $CONDA_ENV_NAME pip install --no-cache-dir --no-build-isolation . \
    && rm -rf /root/.cache/pip \
    && pangolin --version > /pangolin-version.txt
    
WORKDIR /tmp/

ENTRYPOINT ["pangolin"]

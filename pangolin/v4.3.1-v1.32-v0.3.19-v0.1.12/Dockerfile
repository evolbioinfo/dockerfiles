# Largely inspired from
#  https://github.com/cov-lineages/pangolin/blob/master/Dockerfile

FROM continuumio/miniconda3:24.11.1-0
LABEL version="4.3.1" \
      description="Docker image for Pangolin"

ENV VERSION "v4.3.1"
ENV DATAVERSION "v1.32"
ENV SCORPIOVERSION "v0.3.19"
ENV CONSTELLATIONVERSION "v0.1.12"

COPY pangolin.patch /

# Install git for pangolin
RUN apt-get update -y \
    && apt-get install -y git bash libarchive-dev \
    && pip install --upgrade pip \
    && pip uninstall --yes pyOpenSSL \
    && pip install --upgrade cryptography \
    && pip install pyOpenSSL \
    && conda update -y conda \
    && conda update -y --all

RUN cd / \
    && git clone https://github.com/cov-lineages/pangolin.git \
    && cd pangolin \
    && git checkout $VERSION \
    && git apply /pangolin.patch \
    && cd .. \
    && sed -i 's/name: pangolin/name: base/' /pangolin/environment.yml \
    && sed -i s=git+https://github.com/cov-lineages/pangolin-data.git=git+https://github.com/cov-lineages/pangolin-data.git@$DATAVERSION= /pangolin/environment.yml \
    && sed -i s=git+https://github.com/cov-lineages/scorpio.git=git+https://github.com/cov-lineages/scorpio.git@$SCORPIOVERSION= /pangolin/environment.yml \
    && sed -i s=git+https://github.com/cov-lineages/constellations.git=git+https://github.com/cov-lineages/constellations.git@$CONSTELLATIONVERSION= /pangolin/environment.yml \
    && sed -i s/pip=.*/pip/ /pangolin/environment.yml \
    && conda env update --quiet -f /pangolin/environment.yml && conda clean -a \
    && cd /pangolin \
    && pip install . && rm -rf /root/.cache/pip \
    && pangolin --version > /pangolin-version.txt
    
WORKDIR /tmp/

ENTRYPOINT ["pangolin"]

# Largely inspired from
#  https://github.com/cov-lineages/pangolin/blob/master/Dockerfile

FROM continuumio/miniconda3:4.9.2-alpine
LABEL version="4.0.6" \
      description="Docker image for Pangolin"

ENV PATH "/opt/conda/envs/pangolin/bin:$PATH"
ENV VERSION "v4.0.6"
ENV DATAVERSION "v1.8"
ENV SCORPIOVERSION "v0.3.17"
ENV CONSTELLATIONVERSION "v0.1.10"

# Install git for pangolin
RUN apk update && \
    apk add git bash

RUN cd / \
    && git clone https://github.com/cov-lineages/pangolin.git \
    && cd pangolin \
    && git checkout $VERSION \
    && cd .. \
    && sed -i s=git+https://github.com/cov-lineages/pangolin-data.git=git+https://github.com/cov-lineages/pangolin-data.git@$DATAVERSION= /pangolin/environment.yml \
    && sed -i s=git+https://github.com/cov-lineages/scorpio.git=git+https://github.com/cov-lineages/scorpio.git@$SCORPIOVERSION= /pangolin/environment.yml \
    && sed -i s=git+https://github.com/cov-lineages/constellations.git=git+https://github.com/cov-lineages/constellations.git@$CONSTELLATIONVERSION= /pangolin/environment.yml \
    && conda env create --quiet -f /pangolin/environment.yml && conda clean -a \
    && cd /pangolin \
    && pip install . && rm -rf /root/.cache/pip \
    && pangolin --version &> /pangolin-version.txt \
    && conda env export --name pangolin > /pangolin.yml 

WORKDIR /tmp/

ENTRYPOINT ["/opt/conda/envs/pangolin/bin/pangolin"]

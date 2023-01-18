# Largely inspired from
#  https://github.com/cov-lineages/pangolin/blob/master/Dockerfile

FROM continuumio/miniconda3:4.9.2-alpine
LABEL version="4.1.3" \
      description="Docker image for Pangolin"

ENV VERSION "v4.1.3"
ENV DATAVERSION "v1.17"
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
    && sed -i 's/name: pangolin/name: base/' /pangolin/environment.yml \
    && sed -i s=git+https://github.com/cov-lineages/pangolin-data.git=git+https://github.com/cov-lineages/pangolin-data.git@$DATAVERSION= /pangolin/environment.yml \
    && sed -i s=git+https://github.com/cov-lineages/scorpio.git=git+https://github.com/cov-lineages/scorpio.git@$SCORPIOVERSION= /pangolin/environment.yml \
    && sed -i s=git+https://github.com/cov-lineages/constellations.git=git+https://github.com/cov-lineages/constellations.git@$CONSTELLATIONVERSION= /pangolin/environment.yml \
    && conda env update  --quiet -f /pangolin/environment.yml && conda clean -a \
    && cd /pangolin \
    && pip install . && rm -rf /root/.cache/pip \
    && pangolin --version &> /pangolin-version.txt
    
WORKDIR /tmp/

ENTRYPOINT ["pangolin"]

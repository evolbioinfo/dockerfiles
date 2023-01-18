# Largely inspired from https://github.com/cov-lineages/pangolin/blob/master/Dockerfile

FROM continuumio/miniconda3:4.9.2-alpine
LABEL version="3.1.16" \
      description="Docker image for Pangolin"

ENV PATH /opt/conda/envs/pangolin/bin:$PATH
ENV VERSION "v3.1.16"

# Install git for pangolin
RUN apk update && \
    apk add git bash

RUN cd / \
    && git clone https://github.com/cov-lineages/pangolin.git \
    && cd pangolin \
    && git checkout $VERSION \
    && cd .. \
    && grep -v "pandas==1.0.1" pangolin/environment.yml > /environment.yml \
    && sed -i "$(grep -n 'python>=3.7' /environment.yml | cut -f1 -d:)d" /environment.yml  \
    && sed -i "$(grep -n 'pip=' /environment.yml | cut -f1 -d:)d" /environment.yml \
    && conda env create --quiet -f /environment.yml && conda clean -a \
    && cd /pangolin \
    && pip install . && rm -rf /root/.cache/pip \
    && pangolin --version &> /pangolin-version.txt \
    && conda env export --name pangolin > /pangolin.yml

WORKDIR /tmp/

ENTRYPOINT ["/opt/conda/envs/pangolin/bin/pangolin"]

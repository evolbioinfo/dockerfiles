# Catch
# https://github.com/broadinstitute/catch
# v1.5.1


FROM continuumio/miniconda3:4.12.0
LABEL version="1.5.1"
LABEL description="Docker image for Catch"
LABEL maintainer="frederic.lemoine@pasteur.fr"

ENV VERSION "v1.5.1"

RUN apt-get update --fix-missing \
    && apt-get install -y git \
    && cd /usr/local/ \
    && git clone https://github.com/broadinstitute/catch.git \
    && cd catch \
    && git checkout $VERSION \
    && conda install -y mamba  -c conda-forge \
    && mamba install -y numpy==1.23.5 scipy==1.9.3 \
    && pip install -e . \
    && cd .. && rm -rf catch

WORKDIR /tmp/

ENTRYPOINT ["catch"]

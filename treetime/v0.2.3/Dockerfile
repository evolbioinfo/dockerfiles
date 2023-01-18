# TreeTime
# https://github.com/neherlab/treetime
# v0.2.3

# base image: Python 2.7.14
FROM python:2.7.14-slim

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

RUN apt-get update --fix-missing \
    && apt-get install -y wget gcc tk8.6 \
    && pip install numpy scipy biopython matplotlib seaborn \
    && cd /usr/local/ \
    && wget -O v0.2.3.tar.gz https://github.com/neherlab/treetime/archive/v0.2.3.tar.gz \
    && tar -xzvf  v0.2.3.tar.gz \
    && cd treetime-0.2.3 \
    && python setup.py install \
    && cp ancestral_reconstruction.py /usr/local/bin/ \
    && cp temporal_signal.py /usr/local/bin/ \
    && cp timetree_inference.py /usr/local/bin/ \
    && cp examples/*.py /usr/local/bin/ \
    && cd /usr/local/ \
    && rm -rf v0.2.3.tar.gz \
    && rm -rf treetime-0.2.3 \
    && apt-get remove -y wget gcc \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["ancestral_reconstruction.py"]

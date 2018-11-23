# PTP
# https://github.com/zhangjiajie/PTP/
# https://species.h-its.org/ptp/
# master branch commit 4bb2dafdd4a80c68f1aeb29479667639c593ad77

# base image: Python 2.7.14
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

ENV PYTHONPATH="/usr/local/PTP/:${PYTHONPATH}"
ENV PATH="/usr/local/PTP/:${PATH}"

RUN apt-get update --fix-missing \
    && apt-get install -y wget python-setuptools python-numpy python-qt4 unzip \
       	       	       	  python-scipy python-mysqldb python-lxml python-matplotlib \
    && cd /usr/local/ \
    && wget https://github.com/zhangjiajie/PTP/archive/4bb2dafdd4a80c68f1aeb29479667639c593ad77.zip \
    && unzip 4bb2dafdd4a80c68f1aeb29479667639c593ad77.zip \
    && mv PTP-4bb2dafdd4a80c68f1aeb29479667639c593ad77 PTP \
    && rm 4bb2dafdd4a80c68f1aeb29479667639c593ad77.zip \
    && python -c "import matplotlib.font_manager" \
    && apt-get remove -y wget unzip\
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["PTP.py"]

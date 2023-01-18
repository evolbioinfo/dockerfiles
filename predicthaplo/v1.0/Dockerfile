# predicthaplo
# http://bmda.cs.unibas.ch/software.html
# v1.0

# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

RUN apt-get update --fix-missing \
    && apt-get install -y wget gcc g++ make liblapack-dev liblapack3 libopenblas-base libopenblas-dev \
    && cd /usr/local/ \
    && wget -O predicthaplo.tar.gz http://bmda.cs.unibas.ch/software/PredictHaplo-1.0.tgz \
    && tar -xzvf predicthaplo.tar.gz \
    && rm -rf predicthaplo.tar.gz \
    && cd PredictHaplo-1.0 \
    && tar xfvz scythestat-1.0.3.tar.gz \
    && cd scythestat-1.0.3 \
    && ./configure --prefix=/usr/local/PredictHaplo-1.0/NEWSCYTHE \
    && make install \
    && cd .. \
    && make \
    && mv PredictHaplo /usr/local/bin/ \
    && cd /usr/local \
    && rm -rf /usr/local/PredictHaplo-1.0 \
    && apt-get remove -y wget gcc g++ make libopenblas-dev liblapack-dev \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/PredictHaplo"]

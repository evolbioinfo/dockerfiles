# mrbayes
# http://mrbayes.sourceforge.net/
# v3.2.6

# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

RUN apt-get update --fix-missing \
    && apt-get install -y wget gcc g++ make subversion automake autoconf libtool ocl-icd-opencl-dev \
    && cd /usr/local/ \
    && wget https://github.com/beagle-dev/beagle-lib/archive/beagle_release_2_1_2.tar.gz \
    && tar -xzvf beagle_release_2_1_2.tar.gz \
    && cd beagle-lib-beagle_release_2_1_2 \
    && ./autogen.sh && ./configure && make && make install \
    && cd ../ \
    && rm -rf beagle* \
    && svn export https://svn.code.sf.net/p/mrbayes/code/tags/v3.2.6/mrbayes-3.2.6.tar.gz \
    && tar -xzvf mrbayes-3.2.6.tar.gz \
    && cd mrbayes-3.2.6/src \
    && autoconf && ./configure && make && make install \
    && cd /usr/local \
    && rm -rf mrbayes* \
    && apt-get remove -y wget gcc g++ make subversion automake autoconf \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["mb"]

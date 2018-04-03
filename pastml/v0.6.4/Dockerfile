FROM ubuntu:xenial

# Python3 needs an UTF-8 locale, http://bugs.python.org/issue19846
ENV LANG C.UTF-8

# Build-time environmental variable so that apt doesn't complain
ARG DEBIAN_FRONTEND=noninteractive

# Create the development environment
RUN apt-get update --fix-missing \
    && apt-get upgrade -y \
    && apt install -y python3-setuptools python3-dev python3-pip python3 wget \
    && pip3 install --no-cache-dir --upgrade pip

# Copy in the version files
ADD *_VERSION /usr/local/

# File Author / Maintainer
MAINTAINER Anna Zhukova <anna.zhukova@pasteur.fr>

# Install GSL
RUN cd /usr/local/ \
    && wget -O gsl-`cat GSL_VERSION`.tar.gz ftp://ftp.gnu.org/gnu/gsl/gsl-`cat GSL_VERSION`.tar.gz \
    && tar -xzvf gsl-`cat GSL_VERSION`.tar.gz \
    && rm -rf gsl-`cat GSL_VERSION`.tar.gz \
    && cd gsl-`cat GSL_VERSION` \
    && ./configure \
    && make -j 4 \
    && make install \
    && make clean \
    && make distclean \
    && cd /usr/local \
    && rm -rf /usr/local/gsl-`cat GSL_VERSION`

# Install pastml
RUN cd /usr/local/ \
    && pip3 install --no-cache-dir pastml==`cat PASTML_VERSION`
# Install cytopast
RUN cd /usr/local/ \
    && pip3 install --no-cache-dir cytopast==`cat CYTOPAST_VERSION`

# Clean up

# Clean up the development environment to reduce image size
RUN apt remove -y build-essential python3-dev \
    && apt autoremove -y \
    && apt autoclean -y \
    && apt clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

RUN cd /usr/local/ \
    && for vf in $( ls *_VERSION ); do echo "$vf: `cat $vf`"; done

WORKDIR /tmp

# The entrypoint runs cytopast with command line arguments
ENTRYPOINT ["cytopast"]

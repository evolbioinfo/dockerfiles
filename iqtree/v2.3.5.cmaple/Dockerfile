# IQ-TREE
# http://www.iqtree.org/
# 2.3.5 cmaple

# base image: Ubuntu
FROM --platform=linux/amd64 ubuntu:24.04

# File Author / Maintainer
LABEL maintainer="julien.fumey@pasteur.fr"

ENV VERSION="2.3.5.cmaple"
ENV CMAPLEVER="1.0.0"
ENV LSD2="2.4.1"

RUN apt-get update --fix-missing \
    && apt-get install -y wget clang libeigen3-dev libboost-all-dev libomp-dev unzip cmake g++

WORKDIR /usr/local/

RUN wget https://github.com/iqtree/iqtree2/archive/refs/tags/v${VERSION}.zip \
    && unzip v${VERSION}

WORKDIR iqtree2-${VERSION}

RUN cd cmaple/ \
    && wget https://github.com/iqtree/cmaple/archive/refs/tags/v${CMAPLEVER}.tar.gz \
    && tar -xvf v${CMAPLEVER}.tar.gz \
    && mv cmaple-${CMAPLEVER}/* . \
    && rm -rf cmaple-${CMAPLEVER} v${CMAPLEVER}.tar.gz \
    && sed -ie 's/typedef int16_t LengthType;/typedef int32_t LengthType;/g' utils/tools.h \
    && cmake -DCMAKE_C_COMPILER=/usr/bin/gcc -DCMAKE_CXX_COMPILER=/usr/bin/g++ . \
    && make

RUN cd lsd2/ \
    && wget https://github.com/tothuhien/lsd2/archive/refs/tags/v.${LSD2}.zip \
    && unzip v.${LSD2}.zip \
    && mv lsd2-v.${LSD2}/* . \
    && rm -rf lsd2-${LSD2} v.${LSD2}.zip

RUN cmake -DUSE_CMAPLE=ON -DCMAKE_C_COMPILER=/usr/bin/gcc -DCMAKE_CXX_COMPILER=/usr/bin/g++ . \
    && make

RUN mv /usr/local/iqtree2-${VERSION}/iqtree2 /usr/local/bin/iqtree

RUN chmod +x /usr/local/bin/iqtree

WORKDIR /usr/local/

RUN rm -rf v${VERSION}.zip \
    && rm -rf iqtree2-${VERSION} \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur \
    && export PATH=$PATH:/usr/local/bin

RUN ls /usr/local/bin
RUN echo $PATH

ENTRYPOINT ["iqtree"]

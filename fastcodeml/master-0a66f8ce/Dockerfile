# FastCodeML
# https://www2.unil.ch/phylo/bioinformatics/fastcodeml.html
# https://gitlab.isb-sib.ch/phylo/fastcodeml
# master branch, commit 0a66f8cea2786f52e2939555525779d27cbb2a84

# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

RUN apt-get update --fix-missing \
    && apt-get install -y wget git g++ cmake libboost-all-dev libblas-dev liblapack-dev libnlopt-dev mpich \
    && cd /usr/local/ \
    && git clone https://gitlab.isb-sib.ch/phylo/fastcodeml.git \
    && cd fastcodeml \
    && git checkout 0a66f8cea2786f52e2939555525779d27cbb2a84 \
    && sed 's/add_executable(fast ${SRCS})/SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")\nSET(BUILD_SHARED_LIBS OFF)\nSET(CMAKE_EXE_LINKER_FLAGS "-static")\nadd_executable(fast ${SRCS})/' CMakeLists.txt > CMakeLists.tmp \
    && mv CMakeLists.tmp CMakeLists.txt \
    && cmake . \
    && make \
    && cp fast /usr/local/bin \
    && cd ../ \
    && rm -rf /usr/local/fastcodeml \
    && apt-get remove -y wget git g++ cmake libboost-all-dev libblas-dev liblapack-dev libnlopt-dev mpich \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["fast"]

#RAPPAS
#https://github.com/phylo42/RAPPAS
#1.21

# base image: Ubuntu
FROM ubuntu:20.04

# File Author / Maintainer
LABEL maintainer=adrien.le-meur@pasteur.fr

COPY RAPPAS.sh /usr/local/bin/RAPPAS

ENV DEBIAN_FRONTEND="noninteractive" TZ="Europe/Paris"
ENV RAPPASVERSION=1.21
ENV RAPPASCOMMIT=b51eb6a9fea5ccfb9e14388519831dfcff45f901
ENV RAXMLNGVERSION=1.0.2
ENV PHYMLVERSION=v3.3.20200621
ENV PHYMLCOMMIT=2ae8fb0964f12682441df0e78469ca6c0cf1faba
ENV PAMLVERSION=4.8a
ENV PAMLSHORTVERSION=4.8

# Install raxml, phyml, paml then rappas
RUN apt-get update --fix-missing \
	&& apt-get install -y openjdk-8-jdk ant git flex bison libgmp3-dev cmake make g++ wget gcc autoconf automake pkg-config \
	&& update-java-alternatives --set java-1.8.0-openjdk-amd64 \
	&& cd /usr/local/ \
	&& git clone --recursive https://github.com/amkozlov/raxml-ng \
	&& cd raxml-ng \
	&& git checkout ${RAXMLNGVERSION} \
	&& mkdir build && cd build \
	&& cmake .. \
	&& make \
	&& mv ../bin/raxml-ng /usr/local/bin/ \
	&& cd /usr/local/ \
	&& rm -rf raxml-ng \
	&& cd /usr/local/ \
	&& wget -O phyml-${PHYMLVERSION}.tar.gz https://github.com/stephaneguindon/phyml/archive/${PHYMLCOMMIT}.tar.gz \
	&& tar -xzvf phyml-${PHYMLVERSION}.tar.gz \
	&& rm -rf phyml-${PHYMLVERSION}.tar.gz \
	&& cd phyml-${PHYMLCOMMIT} \
	&& sh ./autogen.sh \
	&& ./configure \
	&& make \
	&& make install \
	&& cd /usr/local \
	&& rm -rf /usr/local/phyml-${PHYMLCOMMIT} \
	&& cd /usr/local/ \
	&& wget http://abacus.gene.ucl.ac.uk/software/paml${PAMLVERSION}.macosx.tgz \
	&& tar -xzvf paml${PAMLVERSION}.macosx.tgz \
	&& rm -rf paml${PAMLVERSION}.macosx.tgz \
	&& cd paml${PAMLSHORTVERSION}/src \
	&& make \
	&& cp baseml basemlg chi2 codeml evolver infinitesites mcmctree pamp yn00 /usr/local/bin/ \
	&& cd /usr/local \
	&& rm -rf /usr/local/paml${PAMLSHORTVERSION} \
	&& cd /usr/local/ \
	&& git clone -b master https://github.com/blinard-BIOINFO/RAPPAS.git \
	&& cd RAPPAS \
	&& git checkout b51eb6a9fea5ccfb9e14388519831dfcff45f901 \
	&& ant -f build-cli.xml \
	&& mv /usr/local/RAPPAS/dist/RAPPAS.jar /usr/local/bin/ \
	&& rm -rf /usr/local/RAPPAS \
	&& chmod +x /usr/local/bin/RAPPAS \
	&& apt-get remove -y ant git flex bison libgmp3-dev make cmake g++ wget gcc autoconf automake pkg-config \
	&& apt-get autoremove -y \
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/* \
	&& mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/RAPPAS"]

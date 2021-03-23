#RAPPAS
#https://github.com/phylo42/RAPPAS
#1.21

# base image: Ubuntu
FROM ubuntu:20.04

# File Author / Maintainer
LABEL maintainer=adrien.le-meur@pasteur.fr

COPY RAPPAS.sh /usr/local/bin/RAPPAS

ENV DEBIAN_FRONTEND="noninteractive" TZ="Europe/Paris"

RUN apt-get update --fix-missing \
	&& apt-get install -y openjdk-8-jdk ant git \
	&& update-java-alternatives --set java-1.8.0-openjdk-amd64 \
	&& cd /usr/local/ \
	&& git clone -b master https://github.com/blinard-BIOINFO/RAPPAS.git \
	&& cd RAPPAS \
	&& git checkout b51eb6a9fea5ccfb9e14388519831dfcff45f901 \
	&& ant -f build-cli.xml \
	&& mv /usr/local/RAPPAS/dist/RAPPAS.jar /usr/local/bin/ \
	&& rm -rf /usr/local/RAPPAS \
	&& chmod +x /usr/local/bin/RAPPAS \
	&& apt-get remove -y ant git \
	&& apt-get autoremove -y \
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/* \
	&& mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/RAPPAS"]

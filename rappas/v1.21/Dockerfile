#RAPPAS
#https://github.com/phylo42/RAPPAS
#1.21

# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
LABEL maintainer=adrien.le-meur@pasteur.fr

RUN apt-get update --fix-missing \
	&& apt-get install -y openjdk-8-jdk ant git \
	&& update-java-alternatives --set java-1.8.0-openjdk-amd64 \
	&& cd /usr/local/ \
	&& git clone -b master https://github.com/blinard-BIOINFO/RAPPAS.git \
	&& cd RAPPAS && ant -f build-cli.xml \
	&& mkdir /usr/local/bin/RAPPAS \
	&& mv /usr/local/RAPPAS/dist/RAPPAS.jar /usr/local/bin/RAPPAS \
	&& rm -rf /usr/local/RAPPAS \
	&& printf "#!/bin/bash \n JAVA_OPTS=\"${JAVA_OPTS:='-Xmx2g'}\" \n java $JAVA_OPTS -jar /usr/local/bin/RAPPAS/RAPPAS.jar \${1+\"\$@\"}" > /usr/local/bin/RAPPAS/RAPPAS.sh \
	&& chmod +x /usr/local/bin/RAPPAS/RAPPAS.sh \
	&& apt-get remove -y ant git \
	&& apt-get autoremove -y \
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/* \
	&& mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/RAPPAS/RAPPAS.sh"]

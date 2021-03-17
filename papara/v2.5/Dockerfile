# papara
# https://github.com/sim82/papara_nt
# v2.5

# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
LABEL maintainer=adrien.le-meur@pasteur.fr

RUN apt-get update --fix-missing \
	&& apt-get install -y wget g++ git libboost-all-dev\
	&& cd /usr/local/ \
	&& git clone https://github.com/sim82/papara_nt \
	&& cd papara_nt \
	&& git clone https://github.com/sim82/ivy_mike/ \
	&& sh build_papara2.sh \
	&& cp papara /usr/local/bin/ \
	&& rm -rf /local/usr/papara_nt \
	&& apt-get remove -y wget g++ git \
	&& apt-get autoremove -y \
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/* \
	&& mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/papara"]

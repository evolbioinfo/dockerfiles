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
	&& git checkout 242b717d5e64ff5614dd3f365bcbd809f3d991af \
	&& git clone https://github.com/sim82/ivy_mike/ \
	&& cd ivy_mike/ \
	&& git checkout 3269b7b39dc6c129cfe72708d9086f1e8f8c2c98 \
	&& cd .. \
	&& sh build_papara2.sh \
	&& cp papara /usr/local/bin/ \
	&& rm -rf /local/usr/papara_nt \
	&& apt-get remove -y wget g++ git \
	&& apt-get autoremove -y \
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/* \
	&& mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/papara"]

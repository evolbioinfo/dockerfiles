# MAPLE
# https://github.com/NicolaDM/MAPLE
# v0.7.5

# base image: Ubuntu
FROM pypy:3.10-7.3

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

ENV VERSION="v0.7.5"

RUN git clone https://github.com/NicolaDM/MAPLE.git \
    && cd MAPLE \
    && git checkout MAPLE${VERSION} \
    && cp *.py /usr/local/bin/ \
    && sed -i '1s/^/#\!\/opt\/pypy\/bin\/pypy3\n/' /usr/local/bin/* \
    && chmod +x /usr/local/bin/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/MAPLEv0.7.5.py"]

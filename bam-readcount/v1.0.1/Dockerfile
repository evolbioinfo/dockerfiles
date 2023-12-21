FROM ubuntu:22.04

LABEL maintainer="frederic.lemoine@pasteur.fr"
LABEL maintainer="appliedstochastics@gmail.com"

ENV DEBIAN_FRONTEND noninteractive
ENV DEBCONF_NONINTERACTIVE_SEEN true

ENV VERSION=1.0.1

# Requirements
RUN apt-get update \
    && apt-get install -y build-essential cmake \
    && apt-get install -y git \
    && apt-get clean all \
    && rm -rf /var/lib/apt/lists/* && rm -rf /tmp/* && rm -rf /var/tmp/*

# Clone and build version $VERSION
RUN cd / \
    && git clone https://github.com/genome/bam-readcount \
    && cd bam-readcount \
    && git checkout v$VERSION \
    && rm -rf build \
    && mkdir build \
    && cd build \
    && cmake .. \
    && make

FROM ubuntu:22.04

ENV DEBIAN_FRONTEND noninteractive
ENV DEBCONF_NONINTERACTIVE_SEEN true

LABEL \
    version="v$VERSION" \
    description="Ubuntu 22.04 with bam-readcount v$VERSION"

RUN apt-get update

RUN mkdir -p /opt/bam-readcount/bin
RUN mkdir -p /opt/bam-readcount/test-data
WORKDIR /opt/bam-readcount
COPY --from=0 /bam-readcount/build/bin/bam-readcount /opt/bam-readcount/bin/bam-readcount
COPY --from=0 /bam-readcount/test-data/ref.* /opt/bam-readcount/test-data/
COPY --from=0 /bam-readcount/test-data/test.bam* /opt/bam-readcount/test-data/
RUN ln -s /opt/bam-readcount/bin/bam-readcount /usr/bin/bam-readcount

ENTRYPOINT ["/usr/bin/bam-readcount"]

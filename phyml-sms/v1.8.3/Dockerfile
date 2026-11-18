# PhyML-SMS
# https://gite.lirmm.fr/atgc/sms
# v1.8.3
# The bundled PhyML (20120412) is rebuilt with -fcommon (required by GCC >= 10)

# base image: ubuntu
FROM ubuntu:24.04

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

ENV DEBIAN_FRONTEND="noninteractive" TZ="Europe/Paris"
ENV VERSION=1.8.3
ENV PATH="/usr/local/sms-${VERSION}/:${PATH}"

RUN apt-get update --fix-missing \
    && apt-get install -y --no-install-recommends wget ca-certificates gcc make libc6-dev bc unzip r-base-core \
    && cd /usr/local/ \
    && wget https://gite.lirmm.fr/atgc/sms/-/raw/v${VERSION}/tarball/sms-${VERSION}.zip \
    && unzip sms-${VERSION}.zip \
    && rm -rf sms-${VERSION}.zip \
    && cd sms-${VERSION}/phyml-src \
    && rm -rf phyml phyml-20120412 \
    && tar -xzf phyml-20120412.tar.gz \
    && cd phyml-20120412 \
    && ./configure CFLAGS="-O2 -fcommon" \
    && make \
    && mv src/phyml ../phyml \
    && cd .. \
    && rm -rf phyml-20120412 phyml-20120412.tar.gz \
    && cd /usr/local/sms-${VERSION} \
    && find . -type f -exec chmod 644 {} \; \
    && find . -type d -exec chmod 755 {} \; \
    && find . -name "*.sh" -exec chmod 755 {} \; \
    && find . -name "*.R" -exec chmod 755 {} \; \
    && chmod 755 phyml-src/phyml \
    && apt-get remove -y wget gcc make libc6-dev unzip \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/sms-1.8.3/sms.sh"]

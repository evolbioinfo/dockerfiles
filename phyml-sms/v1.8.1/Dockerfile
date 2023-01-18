# PhyML-SMS
# 
# v1.8.1

# base image: R
FROM r-base:3.4.3

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

ENV PATH="/usr/local/sms-1.8.1/:${PATH}"

RUN apt-get update --fix-missing \
    && apt-get install -y wget gcc make bc \
    && cd /usr/local/ \
    && wget http://www.atgc-montpellier.fr/download/sources/sms/sms-1.8.1.zip \
    && unzip sms-1.8.1.zip \
    && rm -rf sms-1.8.1.zip \
    && cd sms-1.8.1 \
    && make \
    && find . -type f -exec sh -c 'chmod 644 {}' \; \
    && find . -type d -exec sh -c 'chmod 755 {}' \; \
    && find . -name "*.sh" -exec sh -c 'chmod 755 {}' \; \
    && find . -name "*.R" -exec sh -c 'chmod 755 {}' \; \
    && chmod +x phyml-src/phyml \
    && rm -rf phyml-src/phyml-20120412* \
    && cd /usr/local \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/sms-1.8.1/sms.sh"]

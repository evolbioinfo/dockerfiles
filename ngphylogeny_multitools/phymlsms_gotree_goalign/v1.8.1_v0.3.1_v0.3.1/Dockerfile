# PhyML-SMS
# 
# v1.8.1
# goalign
# https://github.com/evolbioinfo/goalign
# v0.3.1
# gotree
# https://github.com/evolbioinfo/gotree
# v0.3.1

# base image: R
FROM evolbioinfo/phyml-sms:v1.8.1

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

ENV PATH="/usr/local/sms-1.8.1/:${PATH}"
ENV PATH=/usr/local/go/bin:/gopath/bin/:$PATH
ENV GOPATH=/gopath

RUN cd /usr/local/bin \
    && wget https://github.com/evolbioinfo/goalign/releases/download/v0.3.1/goalign_amd64_linux \
    && mv goalign_amd64_linux goalign \
    && chmod +x goalign 

RUN cd /usr/local/bin \
    && wget https://github.com/evolbioinfo/gotree/releases/download/v0.3.1/gotree_amd64_linux \
    && mv gotree_amd64_linux gotree \
    && chmod +x gotree

ENTRYPOINT ["/usr/local/sms-1.8.1/sms.sh"]

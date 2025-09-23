# Seqtypedetect
# https://github.com/C3BI-pasteur-fr/ngphylogeny-galaxy/blob/master/tools/seqtypedetect/seqtypedetect.py
# v2.0

# base image: Ubuntu
FROM python:3.13

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

COPY seqtypedetect.py /usr/local/bin/

ENTRYPOINT ["/usr/local/bin/seqtypedetect.py"]

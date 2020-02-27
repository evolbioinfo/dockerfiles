# Sphinx
# http://www.sphinx-doc.org
# v1.8.5

FROM ubuntu:18.04

MAINTAINER Frédéric Lemoine <frederic.lemoine@pasteur.fr>

RUN apt-get -y update \
    && apt-get install -y --no-install-recommends \
        python3 \
        python3-pip \
        python3-setuptools \
        texlive-latex-recommended \
        texlive-fonts-recommended \
        texlive-generic-extra \
        texlive-latex-extra \
	latexmk \
        graphviz\
        python3-graphviz \
	make \
    && python3 -m pip install \
        sphinx==1.8.5 \
        hieroglyph==0.7.1 \
    && apt-get -y autoremove \
    && apt-get -y clean \
    && useradd sphinx

USER sphinx

CMD ["/bin/bash"]

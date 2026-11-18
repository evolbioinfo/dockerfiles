# Sphinx
# http://www.sphinx-doc.org
# v9.2.0
# Sphinx 9.2 requires python >= 3.13: based on the official python image
# (Debian) instead of Ubuntu 18.04. texlive-generic-extra is now
# texlive-plain-generic. hieroglyph is unmaintained (last release 2.1.0 in
# 2020) and may not work with Sphinx 9.

FROM python:3.13-slim

LABEL maintainer="frederic.lemoine@pasteur.fr"

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get -y update \
    && apt-get install -y --no-install-recommends \
        texlive-latex-recommended \
        texlive-fonts-recommended \
        texlive-plain-generic \
        texlive-latex-extra \
        latexmk \
        graphviz \
        make \
    && python3 -m pip install --no-cache-dir \
        sphinx==9.2.0 \
        hieroglyph==2.1.0 \
        graphviz \
    && apt-get -y autoremove \
    && apt-get -y clean \
    && rm -rf /var/lib/apt/lists/* \
    && useradd sphinx \
    && mkdir /pasteur

USER sphinx

CMD ["/bin/bash"]

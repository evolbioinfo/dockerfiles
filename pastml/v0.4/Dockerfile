FROM python:3

# Python3 needs an UTF-8 locale, http://bugs.python.org/issue19846
ENV LANG C.UTF-8

# Copy in the version files
ADD PASTML_VERSION .
ADD CYTOPAST_VERSION .

# Install pastml
RUN pip3 install --no-cache-dir pastml==`cat PASTML_VERSION`
# Install cytopast
RUN pip3 install --no-cache-dir cytopast==`cat CYTOPAST_VERSION`

RUN echo "PASTML version: `cat PASTML_VERSION`"
RUN echo "CYTOPAST version: `cat CYTOPAST_VERSION`"

WORKDIR /tmp

# The entrypoint runs cytopast with command line arguments
ENTRYPOINT ["cytopast"]

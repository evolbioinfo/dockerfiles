FROM ubuntu:xenial

# Python3 needs an UTF-8 locale, http://bugs.python.org/issue19846
ENV LANG C.UTF-8

# Build-time environmental variable so that apt doesn't complain
ARG DEBIAN_FRONTEND=noninteractive

# Create the development environment
RUN apt update && \
    apt install -y python3-setuptools python3-dev python3-pip python3 wget
RUN pip3 install --no-cache-dir --upgrade pip

# Copy in the version files
ADD PASTML_VERSION .
ADD CYTOPAST_VERSION .

# Copy in the main code
ADD cytopast_main.py /
# Install cytopast
RUN pip3 install --no-cache-dir cytopast==`cat CYTOPAST_VERSION`

# Clean up the development environment to reduce image size
RUN apt remove -y build-essential python3.*-dev && \
    apt autoremove -y

# Copy in PASTML and make it executable
RUN mkdir /apps; cd /apps; wget "https://github.com/saishikawa/PASTML/releases/download/`cat ../PASTML_VERSION`/PASTML"
RUN chmod +x /apps/PASTML
ENV PATH="/apps/:${PATH}"

RUN echo "PASTML version: `cat PASTML_VERSION`"
RUN echo "CYTOPAST version: `cat CYTOPAST_VERSION`"

WORKDIR /tmp

# The entrypoint runs cytopast with command line arguments
ENTRYPOINT ["/usr/bin/python3", "/cytopast_main.py"]

FROM python:3.9-slim

RUN mkdir /pasteur

# Install pastml
RUN cd /usr/local/ && pip3 install --no-cache-dir pastml==1.9.42

# The entrypoint runs pastml with command line arguments
ENTRYPOINT ["pastml"]

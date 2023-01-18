FROM evolbioinfo/python:v3.8.2

# Python3 needs an UTF-8 locale, http://bugs.python.org/issue19846
ENV LANG C.UTF-8
ENV LC_ALL C.UTF-8

# Switch to your new user in the docker image
USER root

# Install sierrapy=0.4.1 pandas and biopython
RUN cd /usr/local && pip3 install --no-cache-dir sierrapy==0.4.1 pandas==1.5.2 biopython==1.80

ADD sdrmhunter /usr/local/bin/
RUN chmod +x /usr/local/bin/sdrmhunter

# Switch to your new user in the docker image
USER evolbioinfo

# start working in the "evolbioinfo" home directory
WORKDIR /home/evolbioinfo

# The entrypoint runs sdrmhunter with command line arguments
ENTRYPOINT ["/usr/local/bin/sdrmhunter"]
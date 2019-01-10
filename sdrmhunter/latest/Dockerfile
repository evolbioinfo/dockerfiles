FROM python:3.5-slim

# Python3 needs an UTF-8 locale, http://bugs.python.org/issue19846
ENV LANG C.UTF-8

# Install sierrapy=0.2.1 pandas and biopython
RUN pip3 install --no-cache-dir sierrapy==0.2.1 pandas==0.22.0 biopython==1.70

RUN mkdir /pasteur

ADD sdrm_annotator.py /usr/local/bin/
RUN chmod +x /usr/local/bin/sdrm_annotator.py

# create a non-root user named evolbioinfo,
# give them the password "evolbioinfo" put them in the sudo group
RUN useradd -d /home/evolbioinfo -m -s /bin/bash evolbioinfo && echo "evolbioinfo:evolbioinfo" | chpasswd && adduser evolbioinfo sudo

# start working in the "evolbioinfo" home directory
WORKDIR /home/evolbioinfo

# Make the files owned by evolbioinfo
RUN chown -R evolbioinfo:evolbioinfo /home/evolbioinfo

# Switch to your new user in the docker image
USER evolbioinfo

# The entrypoint runs sdrm_annotator with command line arguments
ENTRYPOINT ["/usr/local/bin/sdrm_annotator.py"]
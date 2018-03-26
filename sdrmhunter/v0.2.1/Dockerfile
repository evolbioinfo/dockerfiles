FROM python:3

# Python3 needs an UTF-8 locale, http://bugs.python.org/issue19846
ENV LANG C.UTF-8

# Install sierrapy=0.2.1 pandas and biopython
RUN pip3 install --no-cache-dir sierrapy==0.2.1 pandas==0.22.0 biopython==1.70

RUN mkdir /pasteur

WORKDIR /app
ADD sdrm_annotator.py .

WORKDIR /tmp

# The entrypoint runs sdrm_annotator with command line arguments
ENTRYPOINT ["python3", "/app/sdrm_annotator.py"]
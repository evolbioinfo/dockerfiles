# python3.6 with packages:
# - ete3==3.1.1
# - biopython==1.72
# - numpy==1.15.4
# - pandas==0.23.4
# - scipy==1.1.0
# - matplotlib==3.0.2
# - hdx-python-country==1.5.4
# - pysam==0.15.2
# - seaborn==0.9.0
# - sierrapy==0.2.1
# - itolapi==3.0.3


FROM evolbioinfo/python:v3.6.9

# Switch to your new user in the docker image
USER root

# Install treetime
RUN cd /usr/local/ && pip3 install --no-cache-dir ete3==3.1.1 biopython==1.72 pandas==0.23.4 scipy==1.1.0 matplotlib==3.0.2 hdx-python-country==1.5.4 numpy==1.15.4 pysam==0.15.2 seaborn==0.9.0 sierrapy==0.2.1 itolapi==3.0.3

# Switch to your new user in the docker image
USER evolbioinfo

# The entrypoint runs pyhton3
ENTRYPOINT ["python3"]

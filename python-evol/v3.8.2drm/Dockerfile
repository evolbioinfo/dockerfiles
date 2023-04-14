# python3.8.2 with packages:
# - ete3==3.1.1
# - biopython==1.76
# - numpy==1.18.4
# - pandas==1.0.3
# - scipy==1.4.1
# - matplotlib==3.2.1
# - hdx-python-country==2.5.5
# - pysam==0.15.4
# - sierrapy==0.4.1
# - wikipedia==1.4.0
# and xz utils

FROM evolbioinfo/python-evol:v3.8.2b

# Switch to your new user in the docker image
USER root

# Install xz-utils

RUN pip3 install --no-cache-dir sierrapy==0.4.1
RUN pip3 install --no-cache-dir wikipedia==1.4.0
RUN pip3 install --no-cache-dir surpyval==0.10.10
RUN pip3 install --no-cache-dir pastml==1.9.40

# Switch to your new user in the docker image
USER evolbioinfo

# The entrypoint runs pyhton3
ENTRYPOINT ["python3"]

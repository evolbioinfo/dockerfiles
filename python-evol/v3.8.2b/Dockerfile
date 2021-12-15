# python3.8.2 with packages:
# - ete3==3.1.1
# - biopython==1.76
# - numpy==1.18.4
# - pandas==1.0.3
# - scipy==1.4.1
# - matplotlib==3.2.1
# - hdx-python-country==2.5.5
# - pysam==0.15.4
# and xz utils

FROM evolbioinfo/python-evol:v3.8.2

# Switch to your new user in the docker image
USER root

# Install xz-utils

RUN apt-get update && apt-get install -y xz-utils

# Switch to your new user in the docker image
USER evolbioinfo

# The entrypoint runs pyhton3
ENTRYPOINT ["python3"]

# python3.7.3 with packages:

FROM python:3.7.3

# Switch to your new user in the docker image
USER root

RUN cd /usr/local/ && pip3 install --no-cache-dir biopython==1.73 ete3==3.1.1 scipy==1.3.1 pandas==0.24.2 numpy==1.16.3 iteround==1.0.2 statsmodels==0.10.1

# The entrypoint runs pyhton3
ENTRYPOINT ["python3"]

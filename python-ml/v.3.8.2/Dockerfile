# python 3.8.2 with 
# - pandas=1.1.1
# - scikit-learn=0.23.2
# - numpy=1.19.1

FROM evolbioinfo/python:v3.8.2

# Switch to root user in the docker image 
USER root

# Install pacakges
RUN cd /usr/local && pip3 install --no-cache-dir \
      pandas==1.1.1 \
      scikit-learn==0.23.2 \
      numpy==1.19.1

# Switch to new user
USER evolbioinfo

# The entrypoint runs python3
ENTRYPOINT ["python3"]

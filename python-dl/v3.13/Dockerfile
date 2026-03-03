FROM python:3.13-slim

# Switch to your new user in the docker image
USER root

# Install phylodeep
RUN cd /usr/local/ && pip3 install torch torchvision --index-url https://download.pytorch.org/whl/cpu
RUN pip3 install pytorch-lightning==2.6.1 tensorflow==2.20 ete3 scipy scikit-learn pandas

# Switch to your new user in the docker image
USER evolbioinfo

# start working in the "evolbioinfo" home directory
WORKDIR /home/evolbioinfo

# The entrypoint runs command line with command line argument
ENTRYPOINT ["python"]
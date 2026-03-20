FROM python:3.12-slim


# Install packages
RUN cd /usr/local/ && pip3 install torch torchvision --index-url https://download.pytorch.org/whl/cpu
RUN pip3 install pytorch-lightning==2.6.1 tensorflow==2.20 ete3 scipy scikit-learn pandas matplotlib seaborn optuna optuna-integration torch-geometric

# Switch to your new user in the docker image
# USER evolbioinfo
#
# # start working in the "evolbioinfo" home directory
# WORKDIR /home/evolbioinfo

# The entrypoint runs command line with command line argument
ENTRYPOINT ["/bin/bash"]
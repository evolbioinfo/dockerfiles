# python3.8.2 with packages:
# - VAMB v3.0.2
# - numpy=1.19.5
# - torch=1.7.1
# - pysam=0.16.0.1
# - torchvision=0.8.2
# - typing_extensions=3.7.4.3
# - Pillow=8.1.0

FROM evolbioinfo/python:v3.8.2

# Switch to your new user in the docker image
USER root

# Install VAMB
RUN apt-get update --fix-missing \
    && apt-get install -y wget gcc \
    && pip3 install numpy==1.19.5 torch==1.7.1 pysam==0.16.0.1 \
    	 torchvision==0.8.2 typing_extensions==3.7.4.3 Pillow==8.1.0 \
	 https://github.com/RasmussenLab/vamb/archive/3.0.2.zip

# Switch to your new user in the docker image
USER evolbioinfo

# The entrypoint runs VAMB
ENTRYPOINT ["vamb"]

# python3.8.2 with packages:
# - VAMB v3.0.2
# - numpy=1.19.5
# - torch=1.7.1
# - pysam=0.16.0.1
# - torchvision=0.8.2
# - typing_extensions=3.7.4.3
# - Pillow=8.1.0

FROM nvidia/cuda:10.2-base-ubuntu18.04

ENV PATH=/miniconda/bin:${PATH}

# Install VAMB
RUN apt-get update --fix-missing \
    && apt-get install -y wget gcc \
    && wget http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && bash Miniconda3-latest-Linux-x86_64.sh -p /miniconda -b \
    && rm Miniconda3-latest-Linux-x86_64.sh \
    && conda update -y conda \
    && conda install -c anaconda -y python=3.8.2 \
    && conda install pip \
    && conda install -c pytorch pytorch torchvision cudatoolkit=10.2 \
    && pip install numpy==1.19.5 torch==1.7.1 pysam==0.16.0.1 \
    	 torchvision==0.8.2 typing_extensions==3.7.4.3 Pillow==8.1.0 \
	 https://github.com/RasmussenLab/vamb/archive/3.0.2.zip

# The entrypoint runs VAMB
ENTRYPOINT ["vamb"]

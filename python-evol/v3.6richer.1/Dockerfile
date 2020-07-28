# python3.6 with packages:
#docutils==0.15.2
#graphql-core==2.3.1
#biopython==1.76
#ete3==3.1.1
#hdx-python-country==2.3.5
#matplotlib==3.1.3
#seaborn==0.10.0
#pastml==1.9.29.4
#phylo-treetime==0.7.4
#pysam==0.15.4
#scipy==1.4.1
#sierrapy==0.2.1
#snakemake==5.10.0
#sympy==1.5.1
#wikipedia==1.4.0
#gql==2.0.0


FROM evolbioinfo/python:v3.6.9

# Switch to your new user in the docker image
USER root

# Make sure commonly needed tools exist
RUN apt-get update --fix-missing \
    && apt-get install -y wget gcc g++ make zip


# Install packages
RUN cd /usr/local/ && pip3 install --no-cache-dir gql==2.0.0 docutils==0.15.2 graphql-core==2.3.1 biopython==1.76 ete3==3.1.1 hdx-python-country==2.3.5 matplotlib==3.1.3 seaborn==0.10.0 pastml==1.9.29.4 phylo-treetime==0.7.4 pysam==0.15.4 scipy==1.4.1 sierrapy==0.2.1 snakemake==5.10.0 sympy==1.5.1 wikipedia==1.4.0

# Switch to your new user in the docker image
USER evolbioinfo

# The entrypoint runs pyhton3
ENTRYPOINT ["python3"]

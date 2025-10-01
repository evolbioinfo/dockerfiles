# Machine Learning Supports
# https://github.com/noaeker/bootstrap_repo
# https://academic.oup.com/bioinformatics/article/40/Supplement_1/i208/7700891
# Commit fec985c
# + raxml-ng v1.2.2

# base image: Ubuntu
FROM python:3.13

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

RUN apt-get install -y git 

RUN pip install numpy pandas ete3 matplotlib networkx biopython scikit-learn

# Install RAxML-NG v1.2.2
ENV RAXMLNGVERSION=1.2.2
RUN apt-get update --fix-missing \
    && apt-get install -y git flex bison libgmp3-dev cmake make g++ libgomp1 procps \
    && cd /usr/local/ \
    && git clone --recursive https://github.com/amkozlov/raxml-ng \
    && cd raxml-ng \
    && git checkout ${RAXMLNGVERSION} \
    && mkdir build && cd build \
    && cmake .. \
    && make \
    && mv ../bin/raxml-ng /usr/local/bin/ \
    && cd /usr/local/ \
    && rm -rf raxml-ng \
    && apt-get remove -y git flex bison libgmp3-dev make cmake g++ \
    && apt-get autoremove -y \
    && apt-get clean

# Install MAD
RUN cd /usr/local/ \
    && wget https://www.mikrobio.uni-kiel.de/de/ag-dagan/ressourcen/mad2-2.zip \
    && unzip mad2-2.zip \
    && cp mad/mad mad/mad.R mad/mad.py /usr/local/bin/ \
    && rm -rf mad mad2-2.zip

# Install Booster
RUN wget -O /usr/local/bin/booster https://github.com/evolbioinfo/booster/releases/download/v0.1.2/booster_linux64 \
    && chmod +x /usr/local/bin/booster

ENV PYTHONPATH=/usr/local/bootstrap_repo

COPY predict_branch_supports /usr/local/bin/
RUN apt-get install -y git \
    && cd /usr/local \
    && git clone https://github.com/noaeker/bootstrap_repo.git \
    && cd bootstrap_repo \
    && git checkout fec985c320b7689bd4d52b067343032889a5977c \
    && sed -i  's/RAXML_NG_EXE.*/RAXML_NG_EXE = \"\/usr\/local\/bin\/raxml-ng\"/g' side_code/config.py \
    && sed -i 's/all_splits.append/all_splits._append/g' simulation_edit/feature_extraction.py \
    && sed -i "s/pattern = '/pattern = r'/g" programs/raxml.py \
    && sed -i "s/re.sub('/re.sub(r'/g" programs/raxml.py \    
    && sed -i "s/pattern = '/pattern = r'/g" side_code/basic_trees_manipulation.py \
    && sed -i 's/pattern = "/pattern = r"/g' side_code/basic_trees_manipulation.py \
    && chmod + /usr/local/bin/predict_branch_supports \
    && apt-get remove -y git

RUN pip install dendropy legacy_cgi

ENTRYPOINT ["/usr/local/bin/predict_branch_supports"]

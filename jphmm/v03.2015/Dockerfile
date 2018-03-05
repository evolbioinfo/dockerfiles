# jpHMM
# http://jphmm.gobics.de/
# v03.2015

# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

ENV JPHMMCIRCOSFILES="'/usr/local/jpHMM/src/circos_input_files/'"
ENV JPHMMPRIORS="'/usr/local/jpHMM/priors/'"

RUN apt-get update --fix-missing \
    && apt-get install -y wget make g++ \
    && cd /usr/local/ \
    && wget -O jpHMM.tar.gz http://jphmm.gobics.de/jpHMM.tar.gz \
    && echo "782764a740bceee073ce719c9eb741e9  jpHMM.tar.gz" | md5sum -c \
    && tar -xzvf jpHMM.tar.gz \
    && rm -f jpHMM.tar.gz \
    && cd jpHMM/src/ \
    && make \
    && ln -s $PWD/determine_UR_and_BPI.pl  $PWD/eval_blat_output.pl  $PWD/extend_circ_gen.pl  $PWD/jphmm_out_to_circos.pl*.pl /usr/local/bin/ \
    && ln -s $PWD/jpHMM /usr/local/bin/ \
    && rm -f ../obj/* \
    && apt-get remove -y wget make g++ \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/jpHMM"]

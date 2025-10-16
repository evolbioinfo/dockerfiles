# Treemmer
# https://github.com/fmenardo/Treemmer
# v0.3

# base image
FROM python:3.13

# File Author / Maintainer
LABEL maintainer="frederic.lemoine@pasteur.fr"

COPY treemmer /usr/local/bin/
RUN chmod +x /usr/local/bin/treemmer

RUN apt-get update && apt-get install -y git

RUN pip3 install joblib ete3 six numpy matplotlib legacy_cgi

WORKDIR /usr/local/

RUN git clone https://github.com/fmenardo/Treemmer.git \
    && cd Treemmer \
    && git checkout a3a1632 \
    && cp -r * /usr/local/bin \
    && cd .. \
    && rm -rf Treemer
  
ENTRYPOINT ["/usr/local/bin/treemmer"]

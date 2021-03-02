# R
# https://cran.r-project.org/
# v4.0.2

# base image: r-base:4.0.2
FROM r-base:4.0.2

# File Author / Maintainer
LABEL maintainer=frederic.lemoine@pasteur.fr

RUN apt-get update && apt-get install -y procps && \
    mkdir /pasteur

ENTRYPOINT ["/usr/bin/R"]

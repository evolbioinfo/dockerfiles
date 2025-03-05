# perl + bioperl + bio-cigar + set-intervaltree + Text::CSV
# https://www.perl.org
# perl v5.32.1

# base image: perl:5.32.1
FROM evolbioinfo/perl:v5.32.1-bioperl

# File Author / Maintainer
LABEL maintainer=frederic.lemoine@pasteur.fr

RUN cpanm Text::CSV

RUN mkdir -p /pasteur

ENTRYPOINT ["/usr/local/bin/perl"]

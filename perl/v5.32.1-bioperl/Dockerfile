# perl + bioperl + bio-cigar + set-intervaltree
# https://www.perl.org
# v5.32.1

# base image: perl:5.32.1
FROM perl:5.32.1

# File Author / Maintainer
LABEL maintainer=frederic.lemoine@pasteur.fr

RUN cpan CJFIELDS/BioPerl-1.7.8.tar.gz \
    && cpan TSIBLEY/Bio-Cigar-1.01.tar.gz \
    && cpan SLOYD/Set-IntervalTree-0.12.tar.gz

RUN mkdir /pasteur

ENTRYPOINT ["/usr/local/bin/perl"]

# MultiQC
# http://multiqc.info/
# v1.35

# base image: official MultiQC image (now published as multiqc/multiqc)
FROM multiqc/multiqc:v1.35

# File Author / Maintainer
LABEL maintainer=frederic.lemoine@pasteur.fr

# The base image runs as the 'multiqc' user
USER root
RUN mkdir /pasteur
USER multiqc

ENTRYPOINT ["multiqc"]

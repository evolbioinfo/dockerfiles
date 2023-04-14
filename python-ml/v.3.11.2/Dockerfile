FROM continuumio/miniconda3:22.11.1

ENV VERSION "v3.11.2"

COPY environment.yml /usr/local/

RUN conda install mamba  -c conda-forge \
    && mamba env update  --quiet -f /usr/local/environment.yml \
    && conda clean -a
    
WORKDIR /tmp/

ENTRYPOINT ["python"]

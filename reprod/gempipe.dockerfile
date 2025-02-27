FROM lazzarigioele/python39_cp:latest

COPY gempipe.yml /home/jovyan/

RUN \
conda config --system --prepend channels bioconda && \
mamba env update -n base --file /tmp/gempipe.yml && \
mamba clean --all -f -y && \
rm /home/jovyan/gempipe.yml




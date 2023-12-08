FROM lazzarigioele/python39:latest

COPY gempipe.yml /tmp/

RUN \
conda config --system --prepend channels bioconda && \
mamba env update -n base --file /tmp/gempipe.yml && \
mamba clean --all -f -y

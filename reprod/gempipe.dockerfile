FROM lazzarigioele/python39_cp:latest

COPY gempipe.yml /tmp/

RUN \
conda config --system --prepend channels bioconda && \
mamba env update -n base --file /tmp/gempipe.yml && \
mamba clean --all -f -y && \
mamba create -n bactabolize -y -c kelwyres -c bioconda -c conda-forge 'bactabolize==1.0.3' 'python=3.9' && \
mamba create -n carveme -y -c bioconda -c conda-forge 'carveme==1.6.2' 'python=3.9'




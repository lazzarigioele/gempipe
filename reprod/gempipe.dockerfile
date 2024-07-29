FROM lazzarigioele/python39_cp:latest

COPY gempipe.yml /tmp/

RUN \
conda config --system --prepend channels bioconda && \
mamba env update -n base --file /tmp/gempipe.yml && \
mamba clean --all -f -y && \
mamba create -n bactabolize -y -c scwatts -c bioconda -c conda-forge 'bactabolize==1.0.2' 




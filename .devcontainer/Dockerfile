FROM --platform=linux/amd64 mambaorg/micromamba:1.5.8 AS base

# Erst mal root werden (falls nicht schon so)
USER root

# Git installieren
RUN apt-get update \
 && apt-get install -y git \
 && rm -rf /var/lib/apt/lists/*

# wieder auf mambauser zurück (wie im micromamba-Image vorgesehen)
USER mambauser

COPY --chown=mambauser:mambauser envs/shiny_env.yml /tmp/envs/shiny_env.yml

RUN ls -R /tmp/envs

RUN micromamba create -y -f /tmp/envs/shiny_env.yml -n shiny && \
    micromamba create -y -n r_packages -c conda-forge -c bioconda \
        r-base r-essentials bioconductor-bayesspace && \
    micromamba clean --all --yes

RUN micromamba run -n r_packages R -q -e ' \
  if (!requireNamespace("devtools", quietly=TRUE)) \
    install.packages("devtools", repos="https://cloud.r-project.org"); \
  devtools::install_github("zhengli09/BASS")'


SHELL ["bash", "-c"]
RUN echo "source /usr/local/bin/_micromamba.sh" >> ~/.bashrc

CMD ["/bin/bash"]

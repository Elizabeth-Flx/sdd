FROM mambaorg/micromamba:1.5.8 AS base

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

RUN micromamba create -y -f /tmp/envs/shiny_env.yml -n shiny  && \
    #micromamba create -y -n bayesspace -c conda-forge r-base r-essentials && \
    #micromamba run -n bayesspace R -q -e ' \
    #  options(repos = c(CRAN = "https://cloud.r-project.org")); \
    #  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager"); \
    #  BiocManager::install("BayesSpace")' && \
    micromamba clean --all --yes

SHELL ["bash", "-c"]
RUN echo "source /usr/local/bin/_micromamba.sh" >> ~/.bashrc

CMD ["/bin/bash"]

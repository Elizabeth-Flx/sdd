FROM --platform=linux/amd64 mambaorg/micromamba:1.5.8 AS base

# Erst mal root werden (falls nicht schon so)
USER root


RUN apt-get update \
 && apt-get install -y \
    git \
    build-essential \
    cmake \
    python3-dev \
 && rm -rf /var/lib/apt/lists/*


# wieder auf mambauser zurück (wie im micromamba-Image vorgesehen)
USER mambauser



# ---- shiny ----
COPY --chown=mambauser:mambauser envs/shiny_env.yml /tmp/envs/shiny_env.yml
RUN micromamba create -y -f /tmp/envs/shiny_env.yml -n shiny


# all tools

# ---- BayesSpace + BASS ----
RUN micromamba create -y -n t_bayesspace_bass -c conda-forge -c bioconda \
    r-base r-essentials bioconductor-bayesspace

RUN micromamba run -n t_bayesspace_bass R -q -e '\
    if (!requireNamespace("devtools", quietly=TRUE)) \
      install.packages("devtools", repos="https://cloud.r-project.org"); \
    devtools::install_github("zhengli09/BASS")'

# ---- IRkernel for Jupyter ----
RUN micromamba run -n t_bayesspace_bass R -q -e '\
    install.packages("IRkernel", repos="https://cloud.r-project.org"); \
    IRkernel::installspec(user = FALSE)'

# ---- GraphST ----
COPY --chown=mambauser:mambauser envs/t_graphst_env.yml /tmp/envs/t_graphst_env.yml
RUN micromamba create -y -f /tmp/envs/t_graphst_env.yml -n t_graphst

# ---- SCAN-IT ----
COPY --chown=mambauser:mambauser envs/t_scanit_env.yml /tmp/envs/t_scanit_env.yml
RUN micromamba create -y -f /tmp/envs/t_scanit_env.yml -n t_scanit

# ---- stLearn ----
COPY --chown=mambauser:mambauser envs/t_stlearn_env.yml /tmp/envs/t_stlearn_env.yml
RUN micromamba create -y -f /tmp/envs/t_stlearn_env.yml -n t_stlearn


RUN micromamba clean --all --yes



SHELL ["bash", "-c"]
RUN echo "source /usr/local/bin/_micromamba.sh" >> ~/.bashrc

CMD ["/bin/bash"]

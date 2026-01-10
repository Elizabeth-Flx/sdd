## 0) CRAN Repo + Timeout
options(repos = c(CRAN = "https://cloud.r-project.org"))
options(timeout = 600000)

## 1) Installer-Tools
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

## 2) SHINY-App Dependencies (u.a. shinyBS)
cran_shiny <- c(
  "shiny", "shinyBS", "shinydashboard", "shinyWidgets",
  "DT", "plotly", "bslib", "htmltools", "htmlwidgets",
  "markdown", "rmarkdown", "ggplot2", "ggrepel", "ggpubr",
  "dplyr", "tidyr", "readr", "stringr", "purrr"
)
install.packages(setdiff(cran_shiny, rownames(installed.packages())))

## 3) EnSDD CRAN Dependencies (aus README + typische Helfer)
cran_ensdd <- c(
  "Seurat", "abind", "mclust", "reticulate", "spdep",
  "igraph", "doParallel", "DR.SC"
)
install.packages(setdiff(cran_ensdd, rownames(installed.packages())))

## 4) Bioconductor Dependencies
bioc_pkgs <- c("SingleCellExperiment", "scran", "scater", "BayesSpace")
BiocManager::install(setdiff(bioc_pkgs, rownames(installed.packages())),
                     ask = FALSE, update = FALSE)

## 5) GitHub Dependencies
remotes::install_github("dmcable/spacexr", build_vignettes = FALSE)

## 6) EnSDD selbst
remotes::install_github("keyalone/EnSDD")

cat("\n✅ Done. Try loading:\n")
library(shinyBS)
library(EnSDD)
cat("✅ shinyBS + EnSDD load OK\n")

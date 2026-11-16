# SysBioSDD: Interactive Exploration and Optimization of Spatial Domain Detection

SysBioSDD is a web platform for running, comparing and combining **spatial domain detection (SDD)** methods on spatial transcriptomics data. Upload an AnnData (`.h5ad`) file, choose which tools to run, and the app gives you:

- the cluster assignment of each tool, plotted on the tissue (optionally over the H&E image);
- an **ensemble consensus clustering** built from all selected tools;
- evaluation metrics (ground truth, spatial coherence, silhouette) and Sankey diagrams comparing the methods.

Every tool runs in its own micromamba environment inside one Docker image, so tools that need incompatible Python or R versions can run side by side.

The full write-up is in [report.pdf](report.pdf).

```
                 ┌─ BASS ───────┐
                 ├─ BayesSpace ─┤
Spatial data ────┼─ GraphST ────┼──► Ensemble ──► Metrics
                 ├─ SCAN-IT ────┤
                 └─ stLearn ────┘
```

## Integrated tools

| Tool | Language | Core method | Uses histology |
|---|---|---|---|
| [BayesSpace](https://github.com/edward130603/BayesSpace) | R | Bayesian model with a Potts spatial prior | no |
| [BASS](https://github.com/zhengli09/BASS) | R | Hierarchical Bayesian model, joint cell type + domain inference | no |
| [stLearn](https://github.com/BiomedicalMachineLearning/stLearn) | Python | GNN on spatial, expression and morphology data | yes |
| [SCAN-IT](https://github.com/zcang/SCAN-IT) | Python | GNN with image-segmentation-style domain extraction | no |
| [GraphST](https://github.com/JinmiaoChenLab/GraphST) | Python | Self-supervised contrastive GNN | no |

BASS runs in every analysis and serves as the reference clustering. The other tools are optional.

## Ensemble method

The ensemble combines the cluster labels of the selected tools into one consensus clustering. It builds on [EnSDD](https://doi.org/10.1038/s42003-024-07001-y) and [SpatialESD](https://doi.org/10.1101/2025.09.11.675735):

1. **Align** all tool outputs to one shared spot index.
2. **Cluster similarity**: compute the Jaccard index between every pair of clusters across all methods.
3. **Random-walk propagation** (t = 10) over the cluster-similarity graph, to capture indirect, transitive relationships between clusters.
4. **Soft spot × spot matrices**: one per method. Spots in the same cluster get maximal similarity, and spots in different but similar clusters get intermediate values (EnSDD uses binary 0/1 instead).
5. **Weighted consensus**: learn one weight per method together with the consensus matrix (weighted Frobenius norm with entropy regularisation on the weights, as in EnSDD).
6. **Final clustering**: hierarchical clustering (average linkage) of the consensus matrix into exactly *K* domains.

These settings were chosen by benchmarking the alternatives. Mutual information instead of Jaccard, Leiden instead of hierarchical clustering, and t = 20 all scored the same or worse (Table 1 of the report). The Leiden variants are still in [ensemble_method/ensemble_helper.py](ensemble_method/ensemble_helper.py).

## Evaluation metrics

| Metric | Type | Better | Needs ground truth |
|---|---|---|---|
| ARI (Adjusted Rand Index) | Agreement with labels | ↑ | yes |
| AMI (Adjusted Mutual Information) | Agreement with labels | ↑ | yes |
| Purity | Agreement with labels | ↑ | yes |
| LISI (Local Inverse Simpson's Index) | Spatial coherence | ↓ | no |
| PAS (Percentage of Abnormal Spots) | Spatial coherence | ↓ | no |
| Silhouette (expression) | Cluster separation in expression space | ↑ | no |
| Silhouette (spatial) | Cluster separation in x/y coordinates | ↑ | no |

## Results

Benchmark: 10x Visium human dorsolateral prefrontal cortex (DLPFC), sample **151671**, with the manually annotated cortical layers as ground truth.

| Method | ARI ↑ | Purity ↑ | AMI ↑ | LISI ↓ | PAS ↓ | Sil. expr ↑ | Sil. spatial ↑ |
|---|---|---|---|---|---|---|---|
| BASS | **0.3864** | **0.5776** | **0.5100** | 0.0918 | 0.0397 | -0.0521 | **0.0835** |
| BayesSpace | 0.2973 | 0.5421 | 0.4210 | 0.1924 | 0.3611 | -0.0458 | -0.0175 |
| GraphST | 0.1572 | 0.4231 | 0.3269 | **0.0237** | **0.0202** | -0.1778 | -0.2084 |
| SCAN-IT | 0.3168 | 0.5277 | 0.4270 | 0.0783 | 0.0253 | **-0.0192** | 0.0612 |
| stLearn | 0.2486 | 0.4601 | 0.3141 | 0.1829 | 0.3628 | -0.0923 | -0.0728 |
| **Ensemble** | 0.3425 | 0.5355 | 0.4804 | 0.0288 | 0.0591 | -0.0882 | -0.0672 |

On ARI and AMI the ensemble ranks second behind BASS. Its spatial coherence (LISI, PAS) is close to the best tool, GraphST. No method reaches strong global separation, so all silhouette scores are close to 0.

For multi-omics data, the ensemble produced the most spatially coherent domains, but its accuracy was lower than that of BayesSpace and stLearn, because finer transcriptional differences were merged. The report discusses this trade-off in section 3.3.

## Multi-omics workflow

RNA-seq and ATAC-seq (SHARE-seq) data can be integrated before spatial domain detection:

```
RNA-seq ─┐                   Spatial reference (MERFISH)
         ├─► SHARE-seq ──► Tangram mapping ──► mapped SHARE-seq ─┐
ATAC-seq ┘                       ▲                               ├─► Seurat integration ──► SDD workflow
                        scRNA-seq┘           ──► mapped scRNA ───┘
```

The Tangram mapping happens before the app is used. In the app you enable **multi-omics** and upload two Tangram-mapped `.h5ad` files (RNA assay and ATAC). [multiomics/run_multiomics_integration.sh](multiomics/run_multiomics_integration.sh) integrates them with Seurat (CCA anchors, `IntegrateData`) and writes a single AnnData file, which then goes through the normal pipeline.

> **Note:** In the current version of `run_multiomics_integration.sh`, the Seurat step is commented out. The script reads a precomputed matrix from `multiomics/integrated_multiomics_matrix2.csv`, which is not in the repository because `*.csv` is in `.gitignore`. To run multi-omics on a fresh clone, re-enable the `seurat.R` call or provide that file yourself.

## Getting started

### Requirements

- Docker (tested with Engine 28.5.1). The image is built for `linux/amd64`.
- Several GB of disk space for the environments, and enough RAM for the tools you select. BayesSpace and BASS in particular are slow on large datasets.

### Option A: VS Code dev container (recommended)

1. Open the repository in VS Code with the *Dev Containers* extension installed.
2. Choose **Reopen in Container**. This builds [.devcontainer/Dockerfile](.devcontainer/Dockerfile), which creates every micromamba environment. The first build takes a while.
3. In the container terminal, start the app:

   ```bash
   micromamba run -n shiny shiny run --launch-browser shiny_app/app.py
   ```

### Option B: plain Docker

The image contains only the environments, not the code, so mount the repository into it:

```bash
docker build -f .devcontainer/Dockerfile -t sysbiosdd .

docker run --rm -it -p 8000:8000 -v "$PWD":/workspace -w /workspace sysbiosdd \
    micromamba run -n shiny shiny run --host 0.0.0.0 --port 8000 shiny_app/app.py
```

Then open <http://localhost:8000>.

### Micromamba environments

| Environment | Defined in | Used for |
|---|---|---|
| `shiny` | [envs/shiny_env.yml](envs/shiny_env.yml) | Web app (Python 3.11, Shiny 1.5.0) |
| `t_bayesspace_bass` | Dockerfile | BayesSpace and BASS (R) |
| `t_graphst` | [envs/t_graphst_env.yml](envs/t_graphst_env.yml) | GraphST |
| `t_scanit` | [envs/t_scanit_env.yml](envs/t_scanit_env.yml) | SCAN-IT |
| `t_stlearn` | [envs/t_stlearn_env.yml](envs/t_stlearn_env.yml) | stLearn, AnnData conversion helpers |
| `sdd_metrics` | [envs/ensemble_metrics.yml](envs/ensemble_metrics.yml) | Ensemble method and metrics |
| `seurat_env` | [envs/seurat_env.yml](envs/seurat_env.yml) | Seurat multi-omics integration |

## Using the app

1. **Load data.** Upload a spatial dataset as `.h5ad`, or tick *Enable multi-omics* and upload the RNA and ATAC files. If `data/adata.h5ad` exists, it is loaded at startup as a demo dataset. Data files are not in the repository, so you have to add it yourself.
2. **Morphology (optional).** If the AnnData object contains histology images, enable *Show morphology (H&E)* and pick an image resolution.
3. **Choose tools.** Select any of BayesSpace, BASS, GraphST, SCAN-IT, stLearn and the ensemble method, and optionally *Calculate metrics after tools run*.
4. **Ground truth (optional).** Select a column of `adata.obs` with known labels to enable ARI, AMI and Purity.
5. **Run** and inspect the per-tool spatial plots, the colourable base spatial plot, the metric plots and the Sankey diagrams. The point size slider helps with dense Visium layouts.

Each session works in its own temporary directory (`$TMPDIR/spatial-shiny/session-<id>/`).

### Input format

- `adata.X`: gene expression counts.
- `adata.obsm["spatial"]`: spot coordinates. Alternatively, `centroid_x`/`centroid_y` or `array_row`/`array_col` in `adata.obs`.
- Optional: histology images in `adata.uns["spatial"]` and a ground-truth label column in `adata.obs`.

## Running steps from the command line

The app calls the scripts below. Run them from the repository root so that the `python -m` module paths resolve.

```bash
# Python tools: AnnData in, one-column cluster CSV out
micromamba run -n t_graphst python -m tools.t_graphst  --input data.h5ad --output csv/graphst.csv --k_clusters 7
micromamba run -n t_scanit  python -m tools.t_scanit   --input data.h5ad --output csv/scanit.csv  --k_clusters 7
micromamba run -n t_stlearn python -m tools.t_stlearn  --input data.h5ad --output csv/stlearn.csv --k_clusters 7

# R tools need the AnnData object exported to CSV first
micromamba run -n t_stlearn python -m tools.adata_to_csv --input data.h5ad --output_dir r_input/
micromamba run -n t_bayesspace_bass Rscript tools/t_bayesspace.R r_input/ csv/bayesspace.csv True 7
micromamba run -n t_bayesspace_bass Rscript tools/t_bass.R       r_input/ csv/bass.csv       True 7

# Ground truth labels (optional, needed for ARI/AMI/Purity)
micromamba run -n t_stlearn python -m tools.dataloader --input data.h5ad --output csv/groundtruth.csv --index_name <obs_column>

# Ensemble: reads every CSV in the directory except groundtruth*/ensemble*, writes csv/ensemble.csv
micromamba run -n sdd_metrics python -m ensemble_method.ensemble_main --input_csvs csv/ --k_clusters 7

# Metrics: writes all_metrics_detailed.csv, all_metrics_summary.csv and plots
micromamba run -n sdd_metrics python -m evaluation_metrics.metrics_main \
    --input_anndata <dir containing the .h5ad> --input_csv csv/ --output metrics/
```

Run `--help` on any of the Python scripts to see all options, for example the clustering method, resolution and graph radius.

## Repository structure

```
├── .devcontainer/       Dockerfile and dev container config (all micromamba envs)
├── envs/                Environment files for each tool
├── tools/               Wrappers for each SDD tool + AnnData/CSV conversion helpers
├── ensemble_method/     Consensus clustering (ensemble_main.py, ensemble_helper.py)
├── evaluation_metrics/  Metric computation and plots (UMAP, line plots, Sankey)
├── multiomics/          RNA + ATAC integration with Seurat
├── shiny_app/           Shiny for Python web app
└── report.pdf           Project report
```

## Limitations and future work

- The benchmark covers a single Visium DLPFC slice. Other slices and technologies (e.g. Xenium) may need different ensemble settings.
- Adding a new tool requires a wrapper script and harmonised output so that it fits the evaluation layer.
- Tools run one after another. Running them in parallel would cut runtime on large datasets.
- Tool hyperparameters are not yet exposed in the interface.
- The ensemble uses only the cluster labels. Using spatial coordinates or morphology directly in the consensus could improve it further.

## Authors

Developed by **Diana Bassow**, **Elizabeth Elodie Lochert**, **Serafina Reif** and **Chiara Thomas** at the Data Science in Systems Biology group, School of Life Sciences, Technical University of Munich. Supervised by Chit Tong Lio and Prof. Markus List.

| Contributor | Contributions |
|---|---|
| Diana Bassow | Shiny app; micromamba Docker container |
| Elizabeth Elodie Lochert | Tool design and architecture; micromamba Docker container |
| Serafina Reif | Ensemble method; metrics and benchmarking |
| Chiara Thomas | Multi-omics integration |

## References

- Zhao et al. *Spatial transcriptomics at subspot resolution with BayesSpace.* Nature Biotechnology (2021).
- Li & Zhou. *BASS: multi-scale and multi-sample analysis enables accurate cell type clustering and spatial domain detection.* Genome Biology (2022).
- Pham et al. *stLearn: integrating spatial location, tissue morphology and gene expression.* bioRxiv (2020).
- Cang et al. *SCAN-IT: Domain segmentation of spatial transcriptomics images by graph neural network.* BMVC (2021).
- Long et al. *Spatially informed clustering, integration, and deconvolution of spatial transcriptomics with GraphST.* Nature Communications (2023).
- Li, Tan & Zhang. *Enhancing spatial domain detection in spatial transcriptomics with EnSDD.* Communications Biology (2024).
- Cao et al. *SpatialESD: spatial ensemble domain detection in spatial transcriptomics.* bioRxiv (2025).
- Maynard et al. *Transcriptome-scale spatial gene expression in the human dorsolateral prefrontal cortex.* Nature Neuroscience (2021).

See [report.pdf](report.pdf) for the full reference list.

library(BayesSpace, warn.conflicts = FALSE, quietly = TRUE)
library(SingleCellExperiment, warn.conflicts = FALSE, quietly = TRUE)
library(Matrix, warn.conflicts = FALSE, quietly = TRUE)

args <- commandArgs(trailingOnly = TRUE)

# Arguments
# args[1] : input directory containing rowData.csv, colData.csv, counts.csv.gz
# args[2] : output file to save clustering results (csv)
# args[3] : do normlog (TRUE/FALSE)
# args[4] : k_clusters (optional, default 12)
# args[5] : pcs for clustering (optional, default 15)

input_dir   <- args[1]
output_file <- args[2]
do_normlog  <- as.logical(args[3])
k_clusters  <- if(length(args) >= 4) as.integer(args[4]) else 12   # default 12
n_pcs       <- if(length(args) >= 5) as.integer(args[5]) else 15  # default 15
gamma       <- if(length(args) >= 6) as.numeric(args[6]) else 2  # default 2


rowData <- read.csv(paste0(input_dir, "/rowData.csv"), stringsAsFactors=FALSE)
colData <- read.csv(paste0(input_dir, "/colData.csv"), stringsAsFactors=FALSE, row.names=1)
counts <-  read.csv(paste0(input_dir, "/counts.csv.gz"),
                   row.names=1, check.names=F, stringsAsFactors=FALSE)

counts <- as.matrix(counts)

sce <- SingleCellExperiment(assays=list(counts=as(counts, "dgCMatrix")),
                            rowData=rowData,
                            colData=colData)

set.seed(102)
sce <- spatialPreprocess(sce, platform="ST", 
                              n.PCs=7, n.HVGs=2000, log.normalize=do_normlog)

colData(sce)$row <- colData(sce)$array_row
colData(sce)$col <- colData(sce)$array_col

set.seed(149)
sce <- spatialCluster(sce, q=k_clusters, platform="ST", d=n_pcs,
                           init.method="mclust", model="t", gamma=gamma,
                           nrep=2000, burn.in=100,
                           save.chain=TRUE)

clusters <- data.frame(
  barcode = colnames(sce),
  spatial_cluster = sce$spatial.cluster
)
write.csv(clusters, output_file, row.names = FALSE)
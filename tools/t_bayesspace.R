library(BayesSpace, warn.conflicts = FALSE, quietly = TRUE)
library(SingleCellExperiment, warn.conflicts = FALSE, quietly = TRUE)
library(Matrix, warn.conflicts = FALSE, quietly = TRUE)

args <- commandArgs(trailingOnly = TRUE)

# Arguments
# args[1] : input directory containing rowData.csv, colData.csv, counts.csv.gz
# args[2] : output file to save clustering results (csv)
# args[3] : do normlog (TRUE/FALSE)

input_dir <- args[1]
output_file <- args[2]
do_normlog <- as.logical(args[3])


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
sce <- spatialCluster(sce, q=12, platform="ST", d=7,
                           init.method="mclust", model="t", gamma=2,
                           nrep=2000, burn.in=100,
                           save.chain=TRUE)

clusters <- data.frame(
  barcode = colnames(sce),
  spatial_cluster = sce$spatial.cluster
)
write.csv(clusters, output_file, row.names = FALSE)
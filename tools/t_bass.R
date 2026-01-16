library(BASS, warn.conflicts = FALSE, quietly = TRUE)

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
storage.mode(counts) <- "numeric"

counts <- counts[, rownames(colData)]

X <- list(counts)

# xy: list of spatial coordinates matrices
xy <- list(
  as.matrix(colData[, c("array_row", "array_col")])
)

BASS_obj <- createBASSObject(
  X = X,
  xy = xy,
  C = 15,          # choose based on biology
  R = 4,           # number of spatial domains
  beta_method = "SW"   # or "fix" for large datasets
)


BASS_obj <- BASS.preprocess(BASS_obj, doLogNormalize = do_normlog, nHVGs = 2000, nPCs = 7)
BASS_obj <- BASS.run(BASS_obj)
BASS_obj <- BASS.postprocess(BASS_obj)

domains <- BASS_obj@results$z[[1]]   # spatial domain labels
colData$spatial_cluster <- domains

colData$barcode <- rownames(colData)
write.csv(colData[c("barcode", "spatial_cluster")], output_file, row.names = FALSE)
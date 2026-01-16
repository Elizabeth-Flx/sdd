library(BASS, warn.conflicts = FALSE, quietly = TRUE)

args <- commandArgs(trailingOnly = TRUE)

# Arguments
# args[1] : input directory containing rowData.csv, colData.csv, counts.csv.gz
# args[2] : output file to save clustering results (csv)
# args[3] : do normlog (TRUE/FALSE)
# args[4] : k_clusters (optional, default 12)
# args[5] : n_topics (optional, default 15)

input_dir   <- args[1]
output_file <- args[2]
do_normlog  <- if(length(args) >= 3) as.logical(args[3]) else TRUE
k_clusters  <- if(length(args) >= 4) as.integer(args[4]) else 12   # default 12
n_topics    <- if(length(args) >= 5) as.integer(args[5]) else 15  # default 15


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
  C = n_topics,          # choose based on biology
  R = k_clusters,           # number of spatial domains
  beta_method = "SW"   # or "fix" for large datasets
)


BASS_obj <- BASS.preprocess(BASS_obj, doLogNormalize = do_normlog)
BASS_obj <- BASS.run(BASS_obj)
BASS_obj <- BASS.postprocess(BASS_obj)

domains <- BASS_obj@results$z[[1]]   # spatial domain labels
colData$spatial_cluster <- domains

colData$barcode <- rownames(colData)
write.csv(colData[c("barcode", "spatial_cluster")], output_file, row.names = FALSE)
library(Seurat)

# Get inputs and out dir via command line arguments
args <- commandArgs(trailingOnly = TRUE)

rna_counts_file <- args[1]
atac_counts_file <- args[2]
output_file <- args[3]

# Load counts from CSV
rna_counts <- read.csv(rna_counts_file, row.names = 1)
atac_counts <- read.csv(atac_counts_file, row.names = 1)

# Ensure common genes
common_genes <- intersect(colnames(rna_counts), colnames(atac_counts))
rna_counts <- rna_counts[, common_genes]
atac_counts <- atac_counts[, common_genes]

# Create Seurat objects for each modality
rna <- CreateSeuratObject(counts = t(rna_counts), assay = "RNA")
atac <- CreateSeuratObject(counts = t(atac_counts), assay = "ATAC")

# Find variable features
rna <- NormalizeData(rna)
rna <- FindVariableFeatures(rna)
atac <- NormalizeData(atac)
atac <- FindVariableFeatures(atac)

# Find integration anchors
anchors <- FindIntegrationAnchors(object.list = list(rna, atac), anchor.features = common_genes)

# Integrate data
integrated <- IntegrateData(anchorset = anchors)

# Extract the integrated expression matrix for RNA cells only
DefaultAssay(integrated) <- "integrated"
integrated_matrix <- as.matrix(GetAssayData(integrated, slot = "data"))

# Save the integrated matrix for downstream tools
write.csv(t(integrated_matrix), file = output_file)
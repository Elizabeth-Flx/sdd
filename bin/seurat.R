library(Seurat)

# --- Updated workflow: Load counts CSVs, create Seurat objects, and output combined gene expression matrix ---

# Load counts from CSV
rna_counts <- read.csv("data/tangram_output/sp_rna_counts.csv", row.names = 1)
atac_counts <- read.csv("data/tangram_output/sp_atac_counts.csv", row.names = 1)

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

# Subset to RNA cells (original RNA cell names)
# rna_cell_names <- colnames(rna)
# common_cells <- intersect(colnames(integrated_matrix), rna_cell_names)
# integrated_matrix_rna <- integrated_matrix[, common_cells, drop = FALSE]

# Save the integrated matrix for downstream tools
write.csv(t(integrated_matrix), file = "data/integrated_rna_matrix.csv")
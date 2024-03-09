# Load the Seurat library
library(Seurat)

# Load the Seurat object
seurat_object <- readRDS(sprintf("../%s/data/%s_discovery.rds", snakemake@params$cellType, snakemake@params$cellType))

# Extract the UMAP data and write to a CSV file
umap_data <- Embeddings(seurat_object, "umap")
write.csv(umap_data, file = snakemake@output$umap)

# Extract the cell state matrix and write to a Matrix Market file
cell_state_matrix <- GetAssayData(seurat_object, slot = "counts")
Matrix::writeMM(cell_state_matrix, file = snakemake@output$cellMat)

# Extract all genes and write to a CSV file
all_genes <- rownames(seurat_object)
write.csv(all_genes, file = snakemake@output$allGenes, row.names = FALSE)

# Extract the seurat clusters
seurat_clusters <- Idents(seurat_object)
write.csv(seurat_clusters, file = snakemake@output$clusters)

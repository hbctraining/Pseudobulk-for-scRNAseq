library(Seurat)
library(tidyverse)

meta <- read.csv("data/meta.csv", row.names=1)

# Celltype IDs have are formatted like: 
# {celltype}_{cluster}
# Removing the underscore
meta$celltype <- sub("_.*", "", meta$cluster_id)
meta <- select(meta, -c(cluster_id))

# The following columns in the metadata have duplicate values:
# - nCount_RNA = nUMI
# - nFeature_RNA = nGene
meta <- select(meta, -c(nUMI, nGene))

# Rename columns for more clarity
meta <- meta %>%
  rename(c("orig.ident"="sample", "sample"="condition"))


# Removing cluster resolutions that will not be used
cols <- c(
  "integrated_snn_res.0.1",
  "integrated_snn_res.0.4",
  "integrated_snn_res.0.6",
  "integrated_snn_res.0.8",
  "integrated_snn_res.1",
  # "integrated_snn_res.1.2",
  "integrated_snn_res.1.4",
  "integrated_snn_res.1.8",
  "seurat_clusters"
)

meta <- meta %>% select(-c(cols))

# Store clusters IDs as factors
meta$seurat_clusters <- meta$integrated_snn_res.1.2
sorted_cluster <- sort(as.integer(unique(meta$seurat_clusters)))
meta$seurat_clusters <- factor(meta$seurat_clusters, levels=sorted_cluster)


library(Seurat)
library(Matrix)
set.seed(1454944673L) # Using the same seed used in the paper

# Load metadata, barcodes, features, and matrix into R
barcodes <- read.csv("data/filtered_counts/barcodes_filtered_raw_counts_for_pseudotime_and_pseudobulk_DGE.tsv", header=FALSE)
features <- read.csv("data/filtered_counts/genes_filtered_raw_counts_for_pseudotime_and_pseudobulk_DGE.tsv", header=FALSE)
counts <- readMM("data/filtered_counts/filtered_raw_counts_for_pseudotime_and_pseudobulk_DGE.mtx")

# Add gene and cell barcode information to count matrix
row.names(counts) <- features$V1
colnames(counts) <- barcodes$V1

# Create seurat object
seurat <- CreateSeuratObject(
  counts, 
  project = "GSE160585", 
  assay = "RNA", 
  meta.data = meta)


# Log normalization
seurat <- NormalizeData(seurat)

# Identify the most variable genes
seurat <- FindVariableFeatures(seurat, 
                               selection.method = "vst",
                               nfeatures = 3000, 
                               verbose = FALSE)


# Split seurat object by sample
split_seurat <- SplitObject(seurat, split.by = "condition")

# Allow R to use more memory
options(future.globals.maxSize = 4000 * 1024^2)

# Run SCTranform on each sample individually
for (i in 1:length(split_seurat)) {
  # Regress out cell cycle scores
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], 
                                   vars.to.regress = c("S.Score", "G2M.Score"), 
                                   vst.flavor = "v2",
                                   variable.features.n = 3000)
}

# Select the most variable features to use for integration
integ_features <- SelectIntegrationFeatures(object.list = split_seurat, 
                                            nfeatures = 3000)
# Prepare the SCT list object for integration
split_seurat <- PrepSCTIntegration(object.list = split_seurat,
                                   anchor.features = integ_features)
# Find best buddies - can take a while to run
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat,
                                        normalization.method = "SCT",
                                        anchor.features = integ_features)
# Integrate across conditions
seurat_integrated <- IntegrateData(anchorset = integ_anchors, 
                                   normalization.method = "SCT")

# Rejoin the layers in the RNA assay that we split earlier
seurat_integrated[["RNA"]] <- JoinLayers(seurat_integrated[["RNA"]])


seurat_integrated <- RunPCA(seurat_integrated, verbose = FALSE)
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:50)
seurat_integrated <- FindNeighbors(seurat_integrated, dims = 1:50)

Idents(seurat_integrated) <- "condition"
DefaultAssay(seurat_integrated) <- "RNA"

saveRDS(seurat_integrated, "data/BAT_GSE160585.rds")
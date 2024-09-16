---
title: "Sample pre-processing"
author: "Noor Sohail"
date: "September 13, 2024"
---

Approximate time: 15 minutes

## Learning Objectives 

* Understand the steps taken to generate the seurat object used for the workshop.

## Sample data

For this workshop, we will be working with a single-cell RNA-seq dataset from [Tseng et al, 2021](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8076094/). The data is available on GEO under the ID [GSE160585](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE160585). The file we will need to create the fully processed seurat object are:

- Metadata csv file
- Counts matrix
- List of features (genes)
- List of cell barcodes

This information comes from a final, completely processed dataset. We have [materials](https://hbctraining.github.io/scRNA-seq_online/) on how to  on how to generate a similarly fully annotated, filtered dataset from single-cell RNA-seq data.

## Pre-processing steps

We have detailed the steps used to generate the seurat object being used for the following lessons. For more details on how to 

1. Download and unzip the dataset from GEO using bash

```bash
#!/bin/bash

# Create data directory to store downloaded files
mkdir -p data/filtered_counts

# Metadata csv file
wget wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE160nnn/GSE160585/suppl/GSE160585%5Fmetadata%5Ffor%5Fpseudotime%5Fand%5Fpseudobulk%5FDGE.csv.gz -O data/meta.csv.gz

# Features, barcodes, and counts matrix
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE160nnn/GSE160585/suppl/GSE160585%5Ffiltered%5Fraw%5Fcounts%5Ffor%5Fpseudotime%5Fand%5Fpseudobulk%5FDGE.tar.gz -O data/filtered_counts.tar.gz

# Unzip and decompress the files
tar -xvf data/filtered_counts.tar.gz -C data/filtered_counts
gunzip data/meta.csv.gz
```


2. Data wrangling of the metadata

```r
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
meta$seurat_clusters <- meta$integrated_snn_res.1.2
```


3. Generate seurat object from files

```r
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
```


4. Log-normalization and highly variable genes

```r
# Log normalization
seurat <- NormalizeData(seurat)

# Identify the most variable genes
seurat <- FindVariableFeatures(seurat, 
                     selection.method = "vst",
                     nfeatures = 3000, 
                     verbose = FALSE)
```


5. SCTransform and regress out cell cycle scores

```r
# Split seurat object by sample
split_seurat <- SplitObject(seurat, split.by = "sample")

# Run SCTranform on each sample individually
for (i in 1:length(split_seurat)) {
    # Regress out cell cycle scores
    split_seurat[[i]] <- SCTransform(split_seurat[[i]], 
                                    vars.to.regress = c("S.Score", "G2M.Score"), 
                                    vst.flavor = "v2",
                                    variable.features.n = 3000)
}
```


6. CCA integration

```r
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
```


7. PCA, nearest neighbors, UMAP

```r
seurat_integrated <- RunPCA(seurat_integrated, verbose = FALSE)
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:50)
seurat_integrated <- FindNeighbors(seurat_integrated, dims = 1:50)
```


8. Save seurat object

```r
saveRDS(seurat_integrated, "data/BAT_GSE160585.rds")
```

---

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

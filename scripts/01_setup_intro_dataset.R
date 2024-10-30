# September 2024
# HBC single-cell RNA-seq DGE workshop

# Single-cell RNA-seq analysis - metadata

# Load libraries
library(Seurat)
library(tidyverse)
library(ggalluvial)

# Load Seurat object
seurat <- readRDS("data/BAT_GSE160585_final.rds")
colnames(seurat@meta.data)


# Number of cells per sample
png("img/sample_info_ncells.png", width=7, height=5, units="in", res=500)
ggplot(seurat@meta.data) +
  geom_bar(aes(x=sample, fill=condition),
           stat="count", color="black") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
  labs(x="Sample name", y="Number of cells")
dev.off()


# UMAPs of condition and sample
png("img/sample_info_sample_umap.png", width=12, height=6, units="in", res=500)
DimPlot(seurat, group.by=c("sample", "condition"))
dev.off()


# UMAP of clusters
png("img/sample_info_cluster_umap.png", width=6, height=6, units="in", res=500)
p <- DimPlot(seurat, group.by="seurat_clusters") + NoLegend()
LabelClusters(p, id = "seurat_clusters",  fontface = "bold", size = 5, bg.colour = "white", bg.r = .2, force = 0)
dev.off()


# Order clusters numerically
order_cluster <- unique(seurat$seurat_clusters) %>% as.numeric() %>% sort() %>% as.character()
seurat$seurat_clusters <- factor(seurat$seurat_clusters, levels=order_cluster)

# Map clusters to celltypes
png("img/sample_info_celltype_map.png", width=15, height=8, units="in", res=500)
ggplot(seurat@meta.data,
       aes(axis1 = seurat_clusters,
           axis2 = celltype,
           fill = celltype)) +
  geom_alluvium() +
  geom_stratum() +
  geom_text(stat = "stratum",
            aes(label=after_stat(stratum))) +
  theme_void() +
  NoLegend() +
  coord_flip()
dev.off()


# UMAP celltype
png("img/sample_info_celltype_umap.png", width=6, height=6, units="in", res=500)
Idents(seurat) <- "celltype"
p <- DimPlot(seurat) + NoLegend()
LabelClusters(p, id = "ident",  fontface = "bold", size = 6,
              bg.colour = "white", bg.r = .2, force = 0)
dev.off()


# Barplot sample proportion by celltype
png("img/sample_info_celltype_prop.png", width=7, height=5, units="in", res=500)
ggplot(seurat@meta.data) +
  geom_bar(aes(x=celltype, fill=condition), 
           position=position_fill(), color="black") +
  theme_classic() +
  labs(x="Celltype", y="Proportion of cells")
dev.off()
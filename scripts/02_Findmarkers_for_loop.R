# September 2024
# HBC single-cell RNA-seq DGE workshop

# Single-cell RNA-seq analysis - FindMarkers

library(Seurat)
library(tidyverse)
library(EnhancedVolcano)


# Load Seurat object
seurat <- readRDS("data/BAT_GSE160585_final.rds")


########## For loop

celltypes <- sort(unique(seurat$celltype))
for (ct in celltypes) {
  
  seurat_ct <- subset(seurat, subset = (celltype == ct))
  # Determine differentiating markers for TN and cold7
  dge_ct <- FindMarkers(seurat_ct,
                        ident.1="cold7",
                        ident.2="TN",
                        test.use="wilcox"
  )
  
  filename <- paste0("results/findmarkers/", ct, ".csv")
  write.csv(dge_ct, filename)
  dge_ct <- dge_ct %>% subset(p_val_adj < 0.05)
  
  
  filename <- paste0("figures/findmarkers/volcano_", ct, ".png")
  png(filename, width=10, height=7, units="in", res=500)
  p <- EnhancedVolcano(dge_ct,
                       row.names(dge_ct),
                       x="avg_log2FC",
                       y="p_val_adj",
                       title=ct)
  print(p)
  dev.off()
  
  
  genes <- rownames(dge_ct)[1:6]
  
  filename <- paste0("figures/findmarkers/genes_vln_", ct, ".png")
  png(filename, width=9, height=6, units="in", res=300)
  p <- VlnPlot(seurat_ct, genes, ncol=3, idents=c("TN", "cold7")) + patchwork::plot_annotation(title = ct)
  print(p)
  dev.off()
  
  
  # Grab the umap coordinates and condition information for each cell
  df <- FetchData(seurat_ct, c("umap_1", "umap_2", "condition"))
  df_tn <- df %>% subset(condition == "TN")
  df_cold7 <- df %>% subset(condition == "cold7")
  
  # Scatterplot of TN cells
  p_tn <- ggplot() +
    geom_point(data=df, aes(x=umap_1, y=umap_2), color="lightgray", alpha=0.5) +
    geom_point(data=df_tn, aes(x=umap_1, y=umap_2), color="#F8766D") +
    theme_classic() +
    ggtitle("TN cells")
  
  # Scatterplot of cold7 cells
  p_cold7 <- ggplot() +
    geom_point(data=df, aes(x=umap_1, y=umap_2), color="lightgray", alpha=0.5) +
    geom_point(data=df_cold7, aes(x=umap_1, y=umap_2), color="#00B8E7") +
    theme_classic() +
    ggtitle("cold7 cells")
  
  # TN and cold7 UMAPs side by side
  filename <- paste0("figures/findmarkers/cells_umap_", ct, ".png")
  png(filename, width=10, height=5, units="in", res=500)
  p <- p_tn + p_cold7 + patchwork::plot_annotation(title = ct)
  print(p)
  dev.off()
  
  
  filename <- paste0("figures/findmarkers/genes_umap_", ct, ".png")
  png(filename, width=9, height=6, units="in", res=300)
  p <- FeaturePlot(seurat_ct, genes, ncol=3) + patchwork::plot_annotation(title = ct)
  print(p)
  dev.off()
  
}
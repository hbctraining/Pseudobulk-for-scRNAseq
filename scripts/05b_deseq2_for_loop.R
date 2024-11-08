########## For loop

celltypes <- sort(unique(seurat$celltype))
for (ct in celltypes) {
  
  # Subset samples by celltype and condition
  bulk_ct <- bulk %>% subset((celltype == ct)  & (condition %in% c("TN", "cold7")))
  
  
  # Number of cells
  filename <- paste0("figures/DESeq2/ncells_", ct, ".png")
  png(filename, width=7, height=5, units="in", res=500)
  p <- ggplot(bulk_ct@meta.data, aes(x=sample, y=n_cells, fill=condition)) +
    geom_bar(stat="identity", color="black") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    labs(x="Sample name", y="Number of cells", title=ct) +
    geom_text(aes(label=n_cells), vjust=-0.5)
  print(p)
  dev.off()
  
  
  # Get count matrix
  cluster_counts <- FetchData(bulk_ct, layer="counts", vars=rownames(bulk_ct))
  # Create DESeq2 object
  # transpose it to get genes as rows
  dds <- DESeqDataSetFromMatrix(t(cluster_counts),
                                colData = bulk_ct@meta.data,
                                design = ~ condition)
  
  
  # Transform counts for data visualization
  filename <- paste0("figures/DESeq2/pca_", ct, ".png")
  png(filename, width=7, height=5, units="in", res=500)
  rld <- rlog(dds, blind=TRUE)
  p <- plotPCA(rld, intgroup=c("condition")) + ggtitle(ct) + theme_classic()
  print(p)
  dev.off()
  
  # Calculate sample correlation
  rld_mat <- assay(rld)
  rld_cor <- cor(rld_mat)
  
  # Change sample names to original values
  # For nicer plots
  rename_samples <- bulk_ct$sample
  colnames(rld_cor) <- str_replace_all(colnames(rld_cor), rename_samples)
  rownames(rld_cor) <- str_replace_all(rownames(rld_cor), rename_samples)
  
  # Plot heatmap
  filename <- paste0("figures/DESeq2/sample_corr_", ct, ".png")
  png(filename, width=7, height=7, units="in", res=500)
  anno <- bulk_ct@meta.data %>%
    select(sample, condition) %>% 
    remove_rownames() %>% 
    column_to_rownames("sample")
  pheatmap(rld_cor, annotation_col=anno, annotation_row=anno, main=ct)
  dev.off()
  
  
  # Run DESeq2 differential expression analysis
  dds <- DESeq(dds)
  plotDispEsts(dds)
  resultsNames(dds)
  contrast <- c("condition", "cold7", "TN")
  
  
  # Results of Wald test
  res <- results(dds, 
                 contrast=contrast,
                 alpha = 0.05)
  
  
  # Shrink the log2 fold changes to be more appropriate using the apeglm method - should cite [paper]() when using this method
  res <- lfcShrink(dds, 
                   coef = "condition_cold7_vs_TN",
                   res=res,
                   type = "apeglm")
  
  filename <- paste0("results/DESeq2/", ct, ".csv")
  write.csv(res, filename)
  
  
  # Set thresholds
  padj.cutoff <- 0.05
  # Turn the results object into a tibble for use with tidyverse functions
  res_tbl <- res %>%
    data.frame() %>%
    rownames_to_column(var="gene") %>% 
    as_tibble()
  
  # Subset the significant results
  sig_res <- dplyr::filter(res_tbl, 
                           padj < padj.cutoff)
  
  
  filename <- paste0("figures/DESeq2/volcano_", ct, ".png")
  png(filename, width=10, height=7, units="in", res=500)
  p <- EnhancedVolcano(sig_res,
                       sig_res$gene,
                       x="log2FoldChange",
                       y="padj",
                       title=ct)
  print(p)
  dev.off()
  
  
  
  # Extract normalized expression for significant genes from the samples
  normalized_counts <- counts(dds, normalized=T) %>% as.data.frame()
  norm_sig <- normalized_counts %>% 
    dplyr::filter(row.names(normalized_counts) %in% sig_res$gene)
  
  # Set a color palette
  heat_colors <- brewer.pal(6, "YlOrRd")
  
  anno <- colData(dds) %>% 
    as.data.frame() %>% 
    select(condition, celltype)
  
  # Run pheatmap using the metadata data frame for the annotation
  filename <- paste0("figures/DESeq2/heatmap_siggenes_", ct, ".png")
  png(filename, width=7, height=10, units="in", res=500)
  pheatmap(norm_sig,
           color = heat_colors,
           cluster_rows = T,
           show_rownames = F,
           annotation = anno,
           border_color = NA,
           fontsize = 10,
           scale = "row", 
           fontsize_row = 10, 
           height = 20,
           main=paste(ct, "significant genes"))
  dev.off()
  
}

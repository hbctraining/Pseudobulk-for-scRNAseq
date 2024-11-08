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

sig_res %>% head()


EnhancedVolcano(sig_res,
                sig_res$gene,
                x="log2FoldChange",
                y="padj"
)


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
pheatmap(norm_sig,
         color = heat_colors,
         cluster_rows = T,
         show_rownames = F,
         annotation = anno,
         border_color = NA,
         fontsize = 10,
         scale = "row", 
         fontsize_row = 10, 
         height = 20)

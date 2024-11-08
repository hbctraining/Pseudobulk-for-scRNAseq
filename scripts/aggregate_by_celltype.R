pb_list <- list()
for (ct in celltypes) {
  
  # Subset cells to one celltype
  seurat_ct <- subset(seurat, subset=(celltype == ct))
  
  # Aggregate to get pseudobulk
  bulk_ct <- AggregateExpression(
              seurat_ct,
              return.seurat = T,
              assays = "RNA",
              group.by = c("celltype", "sample", "condition")
            )
  
  # Add number of cells per sample
  n_cells <- seurat_ct@meta.data %>% 
                dplyr::count(sample, celltype) %>% 
                rename("n"="n_cells")
  n_cells$sample <- str_replace(n_cells$sample, "_", "-")
  meta_bulk_ct <- left_join(bulk_ct@meta.data, n_cells)
  rownames(meta_bulk_ct) <- meta_bulk_ct$orig.ident
  bulk_ct@meta.data <- meta_bulk_ct
  
  pb_list[[ct]] <- bulk_ct
  
}

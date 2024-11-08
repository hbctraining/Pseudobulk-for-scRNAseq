# Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)
plotPCA(rld, intgroup=c("condition")) + theme_classic()

plotPCA(rld, intgroup=c("n_cells")) + theme_classic()


# Calculate sample correlation
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)

# Change sample names to original values
# For nicer plots
rename_samples <- bulk_vsm$sample
colnames(rld_cor) <- str_replace_all(colnames(rld_cor), rename_samples)
rownames(rld_cor) <- str_replace_all(rownames(rld_cor), rename_samples)

# Plot heatmap
anno <- bulk_vsm@meta.data %>%
  select(sample, condition) %>% 
  remove_rownames() %>% 
  column_to_rownames("sample")
pheatmap(rld_cor, annotation_col=anno, annotation_row=anno)

# Run DESeq2 differential expression analysis
dds <- DESeq(dds)

plotDispEsts(dds)
resultsNames(dds)



contrast <- c("condition", "cold7", "TN")

# Results of Wald test
res <- results(dds, 
               contrast=contrast,
               alpha = 0.05)

res %>% head()


# Shrink the log2 fold changes to be more appropriate using the apeglm method - should cite [paper]() when using this method
res <- lfcShrink(dds, 
                 coef = "condition_cold7_vs_TN",
                 res=res,
                 type = "apeglm")


write.csv(res, "results/DESeq2_vsm.csv")
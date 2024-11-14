# Day 2 Answer key

## Functional analysis

**1. Using the code above as a template, run the over-reresentation analysis on the significantly down-regulated genes from the pseudobulk analysis.**

```
# Extract significant results for down-regulated
sigDown <- dplyr::filter(res_tbl_noNAs, padj < 0.05, log2FoldChange < 0)
sigDown_genes <- as.character(sigDown$gene)
# Run GO enrichment analysis 
egoDown <- enrichGO(gene = sigDown_genes, 
                    universe = all_genes,
                    keyType = "SYMBOL",
                    OrgDb = org.Mm.eg.db, 
                    ont = "BP", 
                    pAdjustMethod = "BH", 
                    qvalueCutoff = 0.05, 
                    readable = TRUE)
# Output results from GO analysis to a table
cluster_summaryDown <- data.frame(egoDown)
# Save results
write.csv(cluster_summaryDown, "results/clusterProfiler_VSM_TNvsCold7_downregulated.csv")
```

**How many significant terms do you find?**

```
nrow(cluster_summaryDown)
```

There are 100 GO biological process terms that are downregulated in cold7 vs TN.

**What are some of the prominent biological processes that are observed?**

```
head(cluster_summaryDown$Description)
```

Most of the top terms relate to translation/ribosomes or RNA processing/splicing.

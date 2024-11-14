# Day 2 Answer key

## Functional analysis

**1. Using the code above as a template, run the over-representation analysis on the significantly down-regulated genes from the pseudobulk analysis.**

```r
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

```r
nrow(cluster_summaryDown)
```

There are 100 GO biological process terms that are downregulated in cold7 vs TN.

**What are some of the prominent biological processes that are observed?**

```r
head(cluster_summaryDown$Description)
```

Most of the top terms relate to translation/ribosomes or RNA processing/splicing.

**2. Now that we have run through functional analysis with the results from Pseudobulk DE, let’s see what results we derive from the DGE lists from our FindMarkers DE analysis.**

**Create a significant DE genes data frame from the FM results with an added fold change criteria to reduce the gene list size.** You can do this by running the code below:

```r
sig_fc_dge <- dge_vsm %>% dplyr::filter(p_val_adj < 0.05, abs(avg_log2FC) > 1)
```

**Use this gene list to run over-representation analysis. Be sure to separate genes into up- and down-regulated first. Also keep in mind that the background gene dataset is different than for the DESeq2 analysis.**

```r
# Create background dataset for hypergeometric testing using all tested genes for significance in the results
all_genes_fm <- as.character(rownames(dge_vsm))
# Extract significant results for up- and down-regulated
sigUp_fm <- dplyr::filter(sig_fc_dge, avg_log2FC > 0)
sigUp_fm_genes <- as.character(rownames(sigUp_fm))
sigDown_fm <- dplyr::filter(sig_fc_dge, avg_log2FC < 0)
sigDown_fm_genes <- as.character(rownames(sigDown_fm))
# Run GO enrichment analysis 
egoUp_fm <- enrichGO(gene = sigUp_fm_genes, 
                     universe = all_genes_fm,
                     keyType = "SYMBOL",
                     OrgDb = org.Mm.eg.db, 
                     ont = "BP", 
                     pAdjustMethod = "BH", 
                     qvalueCutoff = 0.05, 
                     readable = TRUE)
egoDown_fm <- enrichGO(gene = sigDown_fm_genes, 
                       universe = all_genes_fm,
                       keyType = "SYMBOL",
                       OrgDb = org.Mm.eg.db, 
                       ont = "BP", 
                       pAdjustMethod = "BH", 
                       qvalueCutoff = 0.05, 
                       readable = TRUE)
# Output results from GO analysis to a table
cluster_summaryUp_fm <- data.frame(egoUp_fm)
cluster_summaryDown_fm <- data.frame(egoDown_fm)
```

**What are the top terms enriched among up-regulated genes?**

```r
head(cluster_summaryUp_fm$Description)
```

We mostly see terms relating to extracellular matrix organization and cell adhesion.

**What are the top terms enriched among down-regulated genes?**

```r
head(cluster_summaryDown_fm$Description)
```

We see terms related to muscle and fat cell development.

**How do these results compare with what we observed from the Pseudobulk DE functional analysis?**


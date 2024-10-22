---
title: "Single-cell RNA-seq: Pseudobulk differential expression analysis"
author: "Noor Sohail, Mary Piper, Lorena Pantano, Amélie Julé, Meeta Mistry, Radhika Khetani"
date: Monday, September 12 2024
---

Approximate time: 40 minutes

## Learning Objectives:

* Sample-level QC
* Understanding dispersion
* Extracting significant results and visualizing it



## Sample quality

A useful initial step in an RNA-seq analysis is to assess overall similarity between samples:

- Which samples are similar to each other, which are different?
- Does this fit the expectation from the experiment’s design?
- What are the major sources of variation in the dataset?

To explore the similarity of our samples, we will be performing sample-level QC using Principal Component Analysis (PCA) and hierarchical clustering methods. Sample-level QC allows us to see how well our replicates cluster together, as well as, observe whether our experimental condition represents the major source of variation in the data. Performing sample-level QC can also identify any sample outliers, which may need to be explored further to determine whether they need to be removed prior to DE analysis.

When using these unsupervised clustering methods, normalization and log2-transformation of the counts improves the distances/clustering for visualization. DESeq2 uses median of ratios method for count normalization and a regularized log transform (rlog) of the normalized counts for sample-level QC as it moderates the variance across the mean, improving the clustering.


### PCA

Principal Component Analysis (PCA) is a dimensionality reduction technique used to emphasize variation and bring out strong patterns in a dataset. Details regarding PCA are given in our [additional materials](https://hbctraining.github.io/DGE_workshop_salmon_online/lessons/03_DGE_QC_analysis.html).

We can run the `rlog()` function from DESeq2 to normalize and rlog transform the raw counts. Then, we can use the `plotPCA()` function to plot the first two principal components. By default, the `plotPCA()` function uses the top 500 most variable genes to compute principal components, but this parameter can be adjusted.

```r
# Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)
plotPCA(rld, intgroup=c("condition")) + theme_classic()
```

<p align="center">
  <img src="../img/pb_pca_condition.png" width="600">
</p>


In this example, we see a nice separation between our samples on PC1 by our condition of interest, which is great; this suggests that our condition of interest is the largest source of variation in our dataset.


It is also useful to check whether the number of cells from which the aggregated counts were derived influences the separation of the samples in the PCA plot. This is particularly useful if you notice an outlier sample, which may be explained by its very low (or very large) cell count compared to others.

```r
plotPCA(rld, intgroup=c("n_cells")) + theme_classic()
```

<p align="center">
  <img src="../img/pb_pca_ncells.png" width="600">
</p>


### Sample correlation

Similar to PCA, hierarchical clustering is another, complementary method for identifying strong patterns in a dataset and potential outliers. The heatmap displays the correlation of gene expression for all pairwise combinations of samples in the dataset. Since the majority of genes are not differentially expressed, samples generally have high correlations with each other (values higher than 0.80). Samples below 0.80 may indicate an outlier in your data and/or sample contamination.

The hierarchical tree can indicate which samples are more similar to each other based on the normalized gene expression values. The color blocks indicate substructure in the data, and you would expect to see your replicates cluster together as a block for each sample group. Additionally, we expect to see samples clustered similar to the groupings observed in a PCA plot.

```r
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
```

<p align="center">
  <img src="../img/pb_sample_corr.png" width="600">
</p>


Now we determine whether we have any outliers that need removing or additional sources of variation that we might want to regress out in our design formula. Since we detected no outliers by PCA or hierarchical clustering, nor do we have any additional sources of variation to regress, we can proceed with running the differential expression analysis.

## Running DESeq2 


Differential expression analysis with DESeq2 involves multiple steps as displayed in the flowchart below in blue. Briefly, DESeq2 will model the **raw counts**, using normalization factors (size factors) to account for differences in library depth. Then, it will estimate the gene-wise dispersions and shrink these estimates to generate more accurate estimates of dispersion to model the counts. Finally, DESeq2 will fit the negative binomial model and perform hypothesis testing using the Wald test or Likelihood Ratio test. All of these steps are explained in detail in our [additional materials](https://hbctraining.github.io/DGE_workshop_salmon_online/schedule/links-to-lessons.html#part-iii-deseq2).

<p align="center">
<img src="../img/de_workflow_salmon_deseq1.png" width="500">
</p>


All of the steps described above are conveniently performed by running the single `DESeq()` function on the DESeq2 object (`dds`) we created earlier.

```r
# Run DESeq2 differential expression analysis
dds <- DESeq(dds)
```

We can check the fit of the model to our data by looking at the plot of dispersion estimates.

```r
plotDispEsts(dds)
```

<p align="center">
<img src="../img/pb_disp_est.png" width="700">
</p>

The plot is encouraging, since we expect our dispersions to decrease with increasing mean and follow the line of best fit.

Now we need to select which comparison we want to make when running DESeq2 by using the `resultsNames()` function.

```r
resultsNames(dds)
```

```
[1] "Intercept"             "condition_cold7_vs_TN"
```

## Results

Now that we have performed the differential expression analysis, we can explore our results for a particular comparison. To denote our comparison of interest, we need to specify the contrasted groups (here, `cold7` vs. `TN). 

Then, we need to perform shrinkage of the log2 fold changes to correct for the fact that the baseline expression level of a gene affects its estimated fold change (for a given gene, a difference in the average read counts between the 2 contrasted groups of 10 will have a greater impact if the baseline expression level of this gene is 20 than if it is 500; therefore, lowly expressed genes are more likely to show inflated log2 fold change values). Here, we use the apeglm method ([Zhu et al., 2018](https://doi.org/10.1093/bioinformatics/bty895)) for shrinkage estimator calculations. Alternative options for shrinkage estimation and the papers to cite if you use them are further described in the [DESeq2 vignette](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#altshrink).

```r
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

res %>% head()
```

```
        baseMean log2FoldChange     lfcSE     pvalue      padj
       <numeric>      <numeric> <numeric>  <numeric> <numeric>
Xkr4     2.80496     -0.0119983  0.258923 0.86291417 0.9487990
Gm1992   0.00000             NA        NA         NA        NA
Rp1      0.00000             NA        NA         NA        NA
Sox17    6.30657      0.0509278  0.256791 0.57599002 0.8147941
Mrpl15 267.61120     -0.0334181  0.147426 0.78522621 0.9185463
Lypla1 152.10108      0.4480539  0.159921 0.00120728 0.0133217
```

This is a great spot to store the results of the comparison

```r
write.csv(res, "results/VSM_cold7_vs_TN.csv")
```

# Significant genes and visualization

Next, we can filter our table for only the significant genes using a p-adjusted threshold of 0.05

```r
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
```

```
  gene   baseMean log2FoldChange lfcSE        pvalue        padj
  <chr>     <dbl>          <dbl> <dbl>         <dbl>       <dbl>
1 Lypla1   152.            0.448 0.160 0.00121       0.0133     
2 Rrs1      68.3          -0.479 0.221 0.00560       0.0442     
3 Prex2      9.82          2.24  0.578 0.00000500    0.000131   
4 Sulf1     90.8           1.39  0.259 0.00000000531 0.000000295
5 Rpl7    3814.           -0.430 0.115 0.0000468     0.000906   
6 Mcm3      19.9           1.20  0.494 0.000682      0.00835    
```

With these results we can use a few different visualization techniques to explore our results:

- Volcano plot of significant genes
- Heatmap of expression for all significant genes
- Scatterplot of normalized expression of top genes

## Volcano plot 

```r
EnhancedVolcano(sig_res,
        sig_res$gene,
        x="log2FoldChange",
        y="padj"
    )
```

<p align="center">
<img src="../img/pb_volcano.png" width="700">
</p>


## Pseudobulk normalized counts heatmap

```r
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
```

<p align="center">
<img src="../img/pb_sig_heatmap_pseudo.png" height="700">
</p>


## Single-cell normalized counts heatmap

TODO fix

```r
# df <- t(FetchData(seurat_vsm, sig_res$gene, assay="RNA", layer="data"))
# anno <- FetchData(seurat_vsm, c("condition", "celltype"))

# pheatmap(df, 
#     color = heat_colors,
#     cluster_cols = F,
#     cluster_rows = T, 
#     show_rownames = F,
#     show_colnames = F,
#     annotation = anno, 
#     border_color = NA, 
#     fontsize = 10, 
#     scale = "row", 
#     fontsize_row = 10, 
#     height = 20)
```

# Top 6 genes

It is important to take a look at some of the top genes that show up. In particular, it is important to evaluate why these genes showed up in the pseudobulked results and contrast it against the gene expression levels at a single-cell level as well.

We may also be interested in determining the total number of significantly upregulated or downregulated genes above a certain fold change threshold (for example log2 fold change (in absolute value) >0.58, which corresponds to a ~50% increase (or ~30% decrease) in gene expression.


```r
genes <- sig_res %>% 
            arrange(padj) %>% 
            subset(abs(log2FoldChange) > 0.6) %>% 
            head()
genes <- genes$gene
genes
```

```
[1] "Rgs5"  "Mt1"   "Emd"   "Nr4a2" "Cwc25" "Cebpb"
```

## Pseudobulked normalized expression scatterplot

Now that we have identified the significant genes, we can plot a scatterplot of the top 6 significant genes. This plot is a good check to make sure that we are interpreting our fold change values correctly, as well.

```r
plot_list <- list()

for (gene in genes) {
    # Save plotcounts to a data frame object
    d <- plotCounts(dds, gene=gene, intgroup="condition", returnData=TRUE)
    d <- d %>% subset(condition %in% c("cold7", "TN"))

    # Plot the normalized counts for each sample
    p <- ggplot(d, aes(x = condition, y = count, color = condition)) + 
            geom_point(position=position_jitter(w = 0.1,h = 0)) +
            theme_bw() +
            ggtitle(gene) +
            theme(plot.title = element_text(hjust = 0.5)) +
            NoLegend()
    plot_list[[gene]] <- p
}

plot_grid(plotlist=plot_list)
```

<p align="center">
<img src="../img/pb_sig_scatter_pseudo.png" height="500">
</p>


## Single-cell normalized expression violin plot

Ultimately, we are evaluating the gene expression at a single-cell level. Therefore it is prudent to go back to the cellular level to see what these same results look like for each individual cell.

```r
DefaultAssay(seurat_vsm) <- "RNA"
Idents(seurat_vsm) <- "condition"
VlnPlot(seurat_vsm, genes, idents=c("cold7", "TN"))
```

<p align="center">
<img src="../img/pb_sig_vln_sc.png" height="500">
</p>

## UMAP

```r
# Grab the umap coordinates and condition information for each cell
df <- FetchData(seurat_vsm, c("umap_1", "umap_2", "condition"))
df_tn <- df %>% subset(condition == "TN")
df_cold7 <- df %>% subset(condition == "cold7")

p_tn <- ggplot() +
        geom_point(data=df, aes(x=umap_1, y=umap_2), color="lightgray", alpha=0.5) +
        geom_point(data=df_tn, aes(x=umap_1, y=umap_2), color="#F8766D") +
        theme_classic() +
        ggtitle("VSM: TN cells")

p_cold7 <- ggplot() +
        geom_point(data=df, aes(x=umap_1, y=umap_2), color="lightgray", alpha=0.5) +
        geom_point(data=df_cold7, aes(x=umap_1, y=umap_2), color="#00B8E7") +
        theme_classic() +
        ggtitle("VSM: cold7 cells")

p_tn + p_cold7
```

<p align="center">
<img src="../img/pb_tn_cold7_umap.png" height="500">
</p>

```r
FeaturePlot(seurat_vsm, genes, ncol=3)
```

<p align="center">
<img src="../img/pb_sig_umap.png" height="500">
</p>

***

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

---
title: "Single-cell RNA-seq: Pseudobulk differential expression analysis"
author: "Noor Sohail, Mary Piper, Lorena Pantano, Amélie Julé, Meeta Mistry, Radhika Khetani"
date: Monday, September 12 2024
---

Approximate time: 40 minutes

## Learning Objectives:

* Recognize the importance of methods for sample-level QC
* Explain the different steps involved in running `DESeq()`


## Sample-level QC

A useful initial step in an RNA-seq analysis is to assess overall similarity between samples:

- Which samples are similar to each other, which are different?
- Does this fit the expectation from the experiment’s design?
- What are the major sources of variation in the dataset?

To explore the similarity of our samples, we will be performing quality checks using Principal Component Analysis (PCA) and a hierarchical clustering approach. 

<p align="center">
  <img src="../img/de_workflow_salmon_qc.png" width="400">
</p>


Sample-level QC allows us to see **how well our replicates cluster together**, as well as, observe whether our experimental condition represents the major source of variation in the data. Performing sample-level QC can also identify any **sample outliers**, which may need to be explored further to determine whether they need to be removed prior to DE analysis.

When using these unsupervised clustering methods, normalization and log2-transformation of the counts improves the distances/clustering for visualization. DESeq2 uses median of ratios method for count normalization and a regularized log transform (rlog) of the normalized counts for sample-level QC as it moderates the variance across the mean, improving the clustering.


### PCA

Principal Component Analysis (PCA) is a dimensionality reduction technique used to emphasize variation and bring out strong patterns in a dataset. Details regarding PCA are given in our [prepared lesson linked here](https://hbctraining.github.io/DGE_workshop_salmon_online/lessons/principal_component_analysis.html
).

We can run the `rlog()` function from DESeq2 to normalize and rlog transform the raw counts. Then, we can use the `plotPCA()` function to plot the first two principal components. By default, the `plotPCA()` function uses the top 500 most variable genes to compute principal components, but this parameter can be adjusted.

```r
# Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)
plotPCA(rld, intgroup=c("condition")) + theme_classic()
```

<p align="center">
  <img src="../img/pb_pca_condition.png" width="600">
</p>


In this example, we see a nice separation between our samples on PC1 by our condition of interest. This suggests that our condition of interest is the largest source of variation in our dataset. There is also a reasonable amount of within group variation for both TN and cold7 samples, with one of the cold7 samples off on its own in the top right quadrant of the plot.

We can check **whether the number of cells** from which the aggregated counts were derived **influences the separation of the samples** in the PCA plot. This is particularly useful if you notice an outlier sample, which may be explained by its very low (or very large) cell count compared to others. Here, the number of cells does not appear to explain the outlier cold7 sample.

```r
plotPCA(rld, intgroup=c("n_cells")) + theme_classic()
```

<p align="center">
  <img src="../img/pb_pca_ncells.png" width="600">
</p>


### Sample correlation

Similar to PCA, hierarchical clustering is another, complementary method for identifying strong patterns in a dataset and potential outliers. The heatmap displays the correlation of gene expression for all pairwise combinations of samples in the dataset. Since the majority of genes are not differentially expressed, samples generally have high correlations with each other (values higher than 0.80). Samples below 0.80 may indicate an outlier in your data and/or sample contamination.

The hierarchical tree can indicate which samples are more similar to each other based on the normalized gene expression values. The **color blocks indicate substructure in the data**, and you would expect to see your replicates cluster together as a block for each sample group. Additionally, we expect to see samples clustered similar to the groupings observed in a PCA plot.

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


Since we detected no outliers by PCA or hierarchical clustering, nor do we have any additional sources of variation to regress, we can proceed with running the differential expression analysis.

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

Everything from **normalization to linear modeling was carried out by the use of a single function**! This function will print out a message for the various steps it performs:

```r

estimating size factors
estimating dispersions
gene-wise dispersion estimates
mean-dispersion relationship
final dispersion estimates
fitting model and testing
```

### Model fitting 
The aggregated count data generated by RNA-seq exhibits overdispersion (variance > mean) and the statistical distribution used to model the counts needs to account for this. As such, DESeq2 uses a negative binomial distribution to model the RNA-seq counts using the equation below:

<p align="center">
<img src="../img/NB_model_formula.png" width="500">
</p>


The two parameters required are the **size factor**, and the gene-wise **dispersion estimate**; both of which can estimated from the observed data. Next, a generalized linear model (GLM) of the NB family is used to fit the data. Modeling is a mathematically formalized way to approximate how the data behaves given a set of parameters.

* **Size factor** is used to make count values from different samples comparable. This is necessary because different samples may have been sequenced to different depths. DESeq2 uses the [median of ratios method](https://hbctraining.github.io/DGE_workshop_salmon_online/lessons/02_DGE_count_normalization.html#deseq2-normalized-counts-median-of-ratios-method) for computing a size factor for each sample.

* **Dispersion** models the **within-group variability by describing how much the variance deviates from the mean**. A dispersion of 1 would indicate that there is no deviance from the mean (i.e mean == variance). A typical RNA-seq dataset, will exhibit some amount of biological variability present across replicates and so we will always have dispersion values less than one.

| | Effect on dispersion |
|:---:|:---:|
| Variance increases | Dispersion increases |
| Mean expression increases | Dispersion decreases |

We can **check the fit of the model** to our data by looking at the plot of gene-wise dispersion estimates and how they compare with mean expression.

```r
plotDispEsts(dds)
```

<p align="center">
<img src="../img/pb_disp_est.png" width="700">
</p>

* In this plot we have dispersion on the y-axis and mean normalized counts on the x-axis. Each **black dot represents a gene and its initial maximum likelihood dispersion estimate (MLE) given the observed data**. Simply looking at the trend of black dots, we observe an inverse relationship between mean and dispersion.
* The **red line** represents a best fit curve to the gene-wise dispersion estimates. The idea behind fitting a curve to the data is that different genes will have different scales of biological variability, but, **across all genes, there will be a distribution of reasonable estimates of dispersion corresponding to a given mean expression level**.
* The **blue dots** represent **gene-wise dispersion estimates shrunken** towards values predicted by the best fit curve. The amount of shrinkage will depend on how far the initial estimate is from the curve, and that distance will vary depending on total number of replicates.
  * If the initial estimate (black dot) is much lower than the fitted curve then values are shrunken up towards the red line.
  * Dispersion estimates that are slightly above the curve are also shrunk toward the curve for better dispersion estimation.
  * Genes with **extremely high dispersion values are not shrunken** (these genes are shown surrounded by **blue circle**). This is due to the likelihood that the gene does not follow the modeling assumptions and has higher variability than others for biological or technical reasons.

Based on the trends observed in this curve, we can say that there is a good fit of the model to the data.

### Setting up contrasts
Now we need to indicate which two sample classes we are interested in comparing, and we do this by specifying **contrasts**. The contrasts are used as input to the DESeq2 `results()` function to extract the desired results. 

> **_NOTE:_ If we run the `results()` function without specifying `contrast` or `name`, it will return the comparison of the last level of the last variable in the design formula over the first level of this variable.** If the order of levels are not specified, they are ordered alphabetically by DESeq2.

We can use the `resultsNames()` function to guide us on exact arguments to provide when extracting our results:

```r
resultsNames(dds)
```

```
[1] "Intercept"             "condition_cold7_vs_TN"
```

To denote our comparison of interest, we need to specify the contrasted groups (here, `cold7` vs. `TN). 

```r
contrast <- c("condition", "cold7", "TN")
```

We use this contrast to extract results:

```r
# Results of Wald test
res <- results(dds, 
            contrast=contrast,
            alpha = 0.05)
```

### Understanding the results

Now let's take a look at **what information is stored** in the results:

```r
res %>% head()
```

**CHANGE THESE VALUES AND PASTE IN THE ACTUAL RESULTS**

```r
        baseMean log2FoldChange     lfcSE     pvalue      padj
       <numeric>      <numeric> <numeric>  <numeric> <numeric>
Xkr4     2.80496     -0.0119983  0.258923 0.86291417 0.9487990
Gm1992   0.00000             NA        NA         NA        NA
Rp1      0.00000             NA        NA         NA        NA
Sox17    6.30657      0.0509278  0.256791 0.57599002 0.8147941
Mrpl15 267.61120     -0.0334181  0.147426 0.78522621 0.9185463
Lypla1 152.10108      0.4480539  0.159921 0.00120728 0.0133217
```

We should have six columns of information reported for each gene (row). We can use the `mcols()` function to extract information on what the values stored in each column represent:

```r
# Get information on each column in results
mcols(res, use.names=T)
```

* `baseMean`: mean of normalized counts for all samples
* `log2FoldChange`: log2 fold change
* `lfcSE`: standard error
* `stat`: Wald statistic
* `pvalue`: Wald test p-value
* `padj`: BH adjusted p-values
 

Although we Then, we need to perform shrinkage of the log2 fold changes to correct for the fact that the baseline expression level of a gene affects its estimated fold change (for a given gene, a difference in the average read counts between the 2 contrasted groups of 10 will have a greater impact if the baseline expression level of this gene is 20 than if it is 500; therefore, lowly expressed genes are more likely to show inflated log2 fold change values). Here, we use the apeglm method ([Zhu et al., 2018](https://doi.org/10.1093/bioinformatics/bty895)) for shrinkage estimator calculations. Alternative options for shrinkage estimation and the papers to cite if you use them are further described in the [DESeq2 vignette](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#altshrink).


```r

# Shrink the log2 fold changes to be more appropriate using the apeglm method - should cite [paper]() when using this method
res <- lfcShrink(dds, 
                coef = "condition_cold7_vs_TN",
                res=res,
                type = "apeglm")
```




**INCLUDE A BRIEF OVERVIEW OF GENE-WISE FILTERING THAT IS APPLIED.**


This is a great spot to store the results of the comparison

```r
write.csv(res, "results/VSM_cold7_vs_TN.csv")
```



***

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

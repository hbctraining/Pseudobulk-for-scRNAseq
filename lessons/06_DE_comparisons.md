## Comparing results from different DE approaches

Approximate time: 40 minutes

## Learning Objectives:

* Compare and contrast results from `DESeq2` and `FindMarkers`
* Evaluate why there are differences in results

## Comparing DGE results


As usual, let's open a new Rscript file, and start with some comments to indicate what this file is going to contain:

```r
# September 2024
# HBC single-cell RNA-seq DGE workshop

# Single-cell RNA-seq analysis - compare DGE results
```

Next we will load the necessary libraries as well.

```r
library(Seurat)
library(tidyverse)
# install.packages("ggvenn")
library("ggvenn")
```

### Load previous results

To start, let us load the results from the the DESeq2 and FindMarkers lessons.

```r
dge_fm <- read.csv("results/findmarkers_vsm.csv")
dge_deseq2 <- read.csv("results/DESeq2_vsm.csv")
```

### Common significant genes

We can visualize how many genes can be found in commmon or are unique to each method by representing the significant genes as a venn diagram.

```r
# Subset to significant genes
sig_fm <- dge_fm %>% subset(p_val_adj < 0.05)
sig_deseq2 <- dge_deseq2 %>% subset(padj < 0.05)

# Create list of just significant gene names 
sig_genes <- list(
  FindMarkers = sig_fm$X,
  DESeq2 = sig_deseq2$X
)

# Create venn diagram
ggvenn(sig_genes)
```

<p align="center">
  <img src="../img/DE_venn.png" width="600">
</p>

For a more conservative approach, we could continue our downstream analysis by using just the conserved significant genes for both conditions. First let us identify which genes those are by first merging and cleaning our dataset.

```r
# Merge FindMarkers and DESeq2 results together
dge <- merge(dge_fm, dge_deseq2, by="X")

# Rename columns to easily understand where results came from
# Remove columns we will not be using
dge <- dge %>% rename("gene"="X") %>%
            rename("padj_fm"="p_val_adj", "padj_deseq2"="padj") %>%
            rename("log2FC_fm"="avg_log2FC", "log2FC_deseq2"="log2FoldChange") %>%
            select(-c("p_val", "baseMean", "lfcSE", "pvalue"))

# Create a column called sig
# Identifies which methods a gene is significant in
dge <- mutate(dge, sig = case_when(
                  ((padj_fm < 0.05) & (padj_deseq2 < 0.05)) ~ "both",
                  (padj_fm < 0.05) ~ "FindMarkers",
                  (padj_deseq2 < 0.05) ~ "DESeq2",
                  ((padj_fm > 0.05) & (padj_deseq2 > 0.05)) ~ "Not Significant"))

dge %>% head()
```

```
           gene log2FC_fm pct.1 pct.2      padj_fm log2FC_deseq2 padj_deseq2
1 0610009B22Rik 0.4919641 0.238 0.138 6.745178e-07   0.178619741   0.5217744
2 0610009O20Rik 0.3421628 0.259 0.152 2.299513e-06   0.021075663   0.9539047
3 0610010K14Rik 0.3505411 0.163 0.088 1.093917e-05  -0.004912021   0.9819494
4 0610012D04Rik 1.1632519 0.033 0.010 3.048507e-03   0.610698014   0.1015138
5 0610012G03Rik 0.1019817 0.497 0.379 1.000000e+00   0.040974487   0.8935888
6 0610030E20Rik 0.6008081 0.147 0.077 6.292494e-06   0.058878725   0.8641972
```

We can again visualize how many genes are not significant and the number significant for each method. You'll notice that we have a few genes that listed `NA`, which is the result of DESeq2 giving an NA for the p-value.

```r
ggplot(dge, aes(x=sig, fill=sig)) +
  geom_bar(stat="count", color="black") +
  theme_classic() + NoLegend() +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
  labs(x="Significant", y="Number of genes") +
  geom_label(vjust=-1,stat="count", aes(label=format(after_stat(count))))
```

<p align="center">
  <img src="../img/DE_ncells.png" width="800">
</p>

With all of this information, we can begin comparing the p-values and average log2-fold changes from each method. To start, we will first remove the genes that are not significant in either method to more clearly see the differences.


```r
# Remove non-significant genes
dge <- dge %>% subset(sig != "Not Significant")

# Compare p-values
ggplot(dge, aes(x=padj_deseq2, y=padj_fm, color=sig)) +
  geom_point() +
  labs(x="DESeq2 Padj", y="FindMarkers Padj", title="Padj values") +
  geom_vline(xintercept = 0.05, color="black", linetype="dashed") +
  geom_hline(yintercept = 0.05, color="black", linetype="dashed") +
  theme_classic()
```

<p align="center">
  <img src="../img/DE_padj.png" width="800">
</p>

```r
# Compare average log2-fold change value
  geom_point() +
  labs(x="DESeq2 LFC", y="FindMarkers LFC", 
        title="Average Log2-fold Change") +
  theme_classic()
```

<p align="center">
  <img src="../img/DE_LFC.png" width="800">
</p>


Next we might ask ourselves, what could be the cause of the differences in the results? If we think back to how we generated the pseudo-bulked results we can consider how the number of cells could effect the final results. There are a couple of different senarios to think of:

* A gene that is highly expressed in very few cells will have a high level of expression at the pseudo-bulk level
* Celltypes where there are few cells will have even fewer cells to aggregate on in the 

Therefore, an important metric to consider is the number or percentage of cells that express the genes we are looking at. We have the columns `pct.1` and `pct.2` which represent the proportion of cells in our dataset that belong to `TN` and `cold7` respectively. So let us consider the data with this additional metric in mind.


```r
pct_1 <- ggplot(dge %>% arrange(pct.1), 
                aes(x=log2FC_deseq2, y=log2FC_fm, color=pct.1)) +
          geom_point() +
          labs(x="DESeq2 LFC", y="FindMarkers LFC", title="Average Log2-fold Change") +
          theme_classic()

pct_2 <- ggplot(dge %>% arrange(pct.2), 
                aes(x=log2FC_deseq2, y=log2FC_fm, color=pct.2)) +
  geom_point() +
  labs(x="DESeq2 LFC", y="FindMarkers LFC", title="Average Log2-fold Change") +
  theme_classic()

pct_1 + pct_2
```

<p align="center">
  <img src="../img/DE_LFC_pct.png" width="800">
</p>




***

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

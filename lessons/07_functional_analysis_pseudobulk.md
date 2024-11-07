---
title: "Single-cell RNA-seq: Pseudobulk visualization"
author: "Noor Sohail, Mary Piper, Lorena Pantano, Amélie Julé, Meeta Mistry, Radhika Khetani"
date: Monday, September 12 2024
---

Approximate time: 40 minutes

## Learning Objectives:

*  Describe the theory of how functional enrichment tools yield statistically enriched functions or interactions
*  Discuss functional analysis using over-representation analysis, and functional class scoring
*  Run clusterProfiler methods on significnat genes from pseudobulk DE analysis

## Functional analysis of differentially expressed genes
In the workshop so far, we have run differential expression analysis using two different approaches: 

1. Using `FindMarkers()` and treating individual cells as replicates
2. Aggregating counts from all cells in sample to run pseudobulk DE

An overview of the comparison between results was provided in the [previous lesson](06_DE_comparisons.md). In this lesson we will take **the results from the pseudobulk DE** and run different **functional analysis** methods to obtain some biological insight. When it comes to functional analysis there are **various analyses** that can be done:

- Determine whether there is enrichment of known biological functions, interactions, or pathways
- Identify genes' involvement in novel pathways or networks by grouping genes together based on similar trends
- Use global changes in gene expression by visualizing all genes being significantly up- or down-regulated in the context of external interaction data

Generally for any differential expression analysis, it is useful to interpret the resulting gene lists using freely available web- and R-based tools.  While tools for functional analysis span a wide variety of techniques, they can loosely be categorized into three main types: over-representation analysis, functional class scoring, and pathway topology [[1](https://github.com/hbctraining/In-depth-NGS-Data-Analysis-Course/raw/master/resources/pathway_tools.pdf)]. 

<p align="center">
<img src="../img/pathway_analysis.png" width="600">
</p>

In this lesson, we will walk you through both an over-representation analysis and gene set enrichment analsysi (GSEA) using an R Bioconductor package called [`clusterProfiler`](https://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html)

## Over-representation analysis
Over-representation analysis (ORA) is used to determine which a priori defined gene sets are more present (over-represented) in a subset of “interesting” genes than what would be expected by chance [(Huang et al, 2009)](https://pmc.ncbi.nlm.nih.gov/articles/PMC2615629/). Most genes in the genome have some pre-existing annotation associated with it which has been compiled through a combination of manual curation and computational algorithms. There are a number of existing databases which define genes using a controlled vocabulary and then categorize genes into groups (gene sets) based on shared function, or involvement in a pathway, or presence in a specific cellular location etc. A very commonly used gene annotataion resource is the [Gene Ontology (GO) database](https://geneontology.org/), and is what we will use in our workflow.

We then use those categorizations to assess enrichment amongst our DE gene results and compare this to the enrichment observed in teh larger univers eof genes, to identify whether or not the enrichment observed is significant.

<p align="center">
<img src="../img/go_proportions.png" width="500">
</p>


### Hypergeometric test
The statistical test that will determine whether something is actually over-represented is the *Hypergeometric test*.

Hypergeometric distribution is a probability distribution that describes the probability of some number of genes (k) being associated with "Functional category 1", for all genes in our gene list (n=1000), compared to drawing some number of genes (K) associated with "Functional category 1" from a population of all of the genes in entire genome (N=13,000) [[2](https://en.wikipedia.org/wiki/Hypergeometric_distribution)].

The calculation of probability of k successes follows the formula:

<p align="center">  
<img src="../img/hypergeo.png" width="200">
</p>

This test will result in an adjusted p-value (after multiple test correction) for each category tested.

### Running ORA with clusterProfiler
Now that we know more about what ORA is doing, let's take our significant genes and see if there are any GO terms over-represented that align with what we expect to be hapepening in VSM cells with change of temperature.

Open up a new R script and let's call it `functional-analysis.R`. The first thing we'll do is load the required libraries:

```r
# Load libraries
library(tidyverse)
library(clusterProfiler)
library(org.Mm.eg.db)
```

Next, we will filter genes to remove any which have NA values in the padj column. These are genes that were not tested and so we do not want to consider them in our background set of genes. Once filtered, we create vectors containing our gene symbols for the background and query set of genes. We will query the up and down-regulated gene sets separately, but note that you can also use the entire significant list as input.

```r
## Create background dataset for hypergeometric testing using all tested genes for significance in the results                 
all_genes <- as.character(res_tbl_noNAs$gene)

## Extract significant results for up-regulated
sigUp <- dplyr::filter(res_tbl_noNAs, padj < 0.05, log2FoldChange > 0)
sigUp_genes <- as.character(sigUp$gene)
```

Finally, we can perform the GO enrichment analysis and save the results:

```r
## Run GO enrichment analysis 
ego <- enrichGO(gene = sigUp_genes, 
                universe = all_genes,
                keyType = "SYMBOL",
                OrgDb = org.Mm.eg.db, 
                ont = "BP", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = TRUE)
```


>**NOTE:** The different organisms with annotation databases available to use with for the `OrgDb` argument can be found [here](../img/orgdb_annotation_databases.png).
>
> Also, the `keyType` argument may be coded as `keytype` in different versions of clusterProfiler.
>
> Finally, the `ont` argument can accept either "BP" (Biological Process), "MF" (Molecular Function), and "CC" (Cellular Component) subontologies, or "ALL" for all three.

```r
## Output results from GO analysis to a table
cluster_summary <- data.frame(ego)

write.csv(cluster_summary, "results/clusterProfiler_VSM_TNv.csv")
```          

> **NOTE:** Instead of saving just the results summary from the `ego` object, it might also be beneficial to save the object itself. The `save()` function enables you to save it as a `.rda` file, e.g. `save(ego, file="results/ego.rda")`. The statistsics stored in the object can be used for downstream visualization.


"In VSMs, extracellular matrix organization, angiogenesis, cell division, cell junction assembly, epithelial cell migration, and response to TGFB stimulus were significantly over-represented in the transcripts whose expression was upregulated by cold (Extended Data Figure 3a-b). In the transcripts downregulated in cold, regulation of RNA splicing, ribosome biogenesis, striated muscle development, and G1/S transition of cell cycle were significantly over-represented (Extended Data Figure 3c-d)."



***

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

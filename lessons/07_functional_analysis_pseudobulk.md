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
In the workshop so far, we have run differential expression analysis using two differnt approaches: 

1. Using `FindMarkers()` and treating individual cells as replicates
2. Aggregating counts from all cells in sample to run pseudobulk DE

An overview of the comparison between results was provided in the [previous lesson](06_DE_comparisons.md). In this lesson we will take the results from the pseudobulk DE and run different functional analysis methods to obtain some biological insight.


- use the DE genes from Pseudobulk analysis to perform functional analysis (clusterProfiler)
  - ORA
  - GSEA


***

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

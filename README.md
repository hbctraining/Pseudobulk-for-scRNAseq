# Pseudobulk and related approaches for scRNA-seq analysis

| Audience | Computational skills required| Duration |
:----------|:----------|:----------|
| Biologists | [Introduction to R](https://hbctraining.github.io/Intro-to-R-flipped/) | 3-session online workshop (~7.5 hours of trainer-led time)|

### Description

This repository has teaching materials for a hands-on workshop on **Pseudobulk and related approaches for scRNA-seq analysis**. 
After clustering cells as part of a single-cell RNA-seq experiment, investigators are often interested in carrying out a differential expression analysis between conditions within certain cell types. When the experimental design includes replicates with multiple conditions there is power in identifying differentially expressed genes using a pseudobulk approach. This workshop will:

- Discuss various approaches for differential expression using single cell RNA-seq data
- Prepare single-cell RNA-seq raw count data for pseudobulk analysis
- Perform differential expression analysis on pseudobulk counts using DESeq2
- Present approaches for evaluating differential proportions of cells between conditions

Working knowledge of R is required or completion of the [Introduction to R workshop](https://hbctraining.github.io/Intro-to-R/). 

**Note for Trainers:** Please note that the schedule linked below assumes that learners will spend between 3-4 hours on reading through, and completing exercises from selected lessons between classes. The online component of the workshop focuses on more exercises and discussion/Q & A.

> These materials were developed for a trainer-led workshop, but are also amenable to self-guided learning.

### Learning Objectives

* Understanding considerations for when to use different DGE algorithms on scRNA-seq data
* Using FindMarkers to evaluate significantly DE genes
* Aggregating single cell expression data into a pseudobulk counts matrix to run a DESeq2 workflow
* Evaluating expression patterns of differentially expressed genes at the pseudobulk and single cell level
* Application of methods for evaluating differential proportions of cells between conditions


### Lessons
* [Workshop schedule (trainer-led learning)](schedule/)
* [Self-learning](schedule/self-learning.md)

### Installation Requirements

#### Applications
Download the most recent versions of R and RStudio for your laptop:

 - [R](http://lib.stat.cmu.edu/R/CRAN/) **(version 4.0.0 or above)**
 - [RStudio](https://www.rstudio.com/products/rstudio/download/#download)

#### Packages for R

> **Note 1: Install the packages in the order listed below.**

> **Note 2:  All the package names listed below are case sensitive!**
 
> **Note 3**: **If you have a Mac,** download and install this tool before intalling your packages:
https://mac.r-project.org/tools/gfortran-12.2-universal.pkg

> **Note 4**: At any point (especially if you’ve used R/Bioconductor in the past), in the console **R may ask you if you want to update any old packages by asking Update all/some/none? [a/s/n]:**. If you see this, **type "a" at the prompt and hit Enter** to update any old packages. _Updating packages can sometimes take quite a bit of time to run, so please account for that before you start with these installations._  

> **Note 5:** If you see a message in your console along the lines of “binary version available but the source version is later”, followed by a question, **“Do you want to install from sources the package which needs compilation? y/n”, type n for no, and hit enter**.


**(1)** Install the 4 packages listed below from **Bioconductor** using the the `BiocManager::install()` function.

```
BiocManager::install("DESeq2")
BiocManager::install("EnhancedVolcano")
BiocManager::install("SingleCellExperiment")
BiocManager::install("miloR")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Mm.eg.db")
```

**Please install them one-by-one as follows:**

```r
BiocManager::install("DESeq2")
BiocManager::install("EnhancedVolcano")
& so on ...
```

**(2)** Install the 7 packages listed below from **CRAN** using the `install.packages()` function. 

```
install.packages("Seurat")
install.packages("tidyverse")
install.packages("pheatmap")
install.packages("RColorBrewer")
install.packages("cowplot")
install.packages("dplyr")
install.packages("ggalluvial")
install.packages("msigdbr")
install.packages("ggvenn")
```

**Please install them one-by-one as follows:**

```r
install.packages("Seurat")
install.packages("tidyverse")
install.packages("pheatmap")
& so on ...
```

**(3)** Finally, please check that all the packages were installed successfully by **loading them one at a time** using the `library()` function.  

```r
library(Seurat)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(cowplot)
library(dplyr)
library(DESeq2)
library(EnhancedVolcano)
library(SingleCellExperiment)
library(miloR)
library(ggalluvial)
library(clusterProfiler)
library(org.Mm.eg.db)
library(msigdbr)
library(ggvenn)
```

**(4)** Once all packages have been loaded, run sessionInfo().  

```r
sessionInfo()
```

---

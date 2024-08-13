# Differential expression analysis of scRNAseq

| Audience | Computational skills required| Duration |
:----------|:----------|:----------|
| Biologists | [Introduction to R](https://hbctraining.github.io/Intro-to-R-flipped/) | 3-session online workshop (~7.5 hours of trainer-led time)|

### Description

This repository has teaching materials for a hands-on **Differential expression analysis of single-cell RNA-seq analysis** workshop. This workshop will inform participants on various approaches for differential expression analysis of single cell RNA-seq and provide some discussion on differential abundance of cells. This will be a hands-on workshop in which we will focus on using .....

Working knowledge of R is required or completion of the [Introduction to R workshop](https://hbctraining.github.io/Intro-to-R/). 

**Note for Trainers:** Please note that the schedule linked below assumes that learners will spend between 3-4 hours on reading through, and completing exercises from selected lessons between classes. The online component of the workshop focuses on more exercises and discussion/Q & A.

> These materials were developed for a trainer-led workshop, but are also amenable to self-guided learning.

### Learning Objectives
* LO1
* 

### Lessons
* [Workshop schedule (trainer-led learning)](schedule/)
* [Self-learning](schedule/links-to-lessons.md)

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


**Please install them one-by-one as follows:**

```r
BiocManager::install("AnnotationHub")
BiocManager::install("ensembldb")
& so on ...
```

**(2)** Install the 8 packages listed below from **CRAN** using the `install.packages()` function. 


**Please install them one-by-one as follows:**

```r
install.packages("tidyverse")
install.packages("Matrix")
install.packages("RCurl")
& so on ...
```

**(3)** Finally, please check that all the packages were installed successfully by **loading them one at a time** using the `library()` function.  

```r
library(Seurat)

```

**(4)** Once all packages have been loaded, run sessionInfo().  

```r
sessionInfo()
```

---


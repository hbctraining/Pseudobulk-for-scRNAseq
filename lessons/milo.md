---
title: "Single-cell RNA-seq: MiloR differential abundance analysis"
author: "Noor Sohail, Meeta Mistry"
date: Monday, September 12 2024
---

Approximate time: 90 minutes

### Learning Objectives:

* Utilize the MiloR tool to perform  differential abundance analysis on a specific cell type cluster


* * *

# Differential abundance analysis with MiloR

Looking at single-cell datasets on a cluster/celltype level is a very common mode of analysis. However, perhaps you have questions on the more subtle shifts within a certain cell population. The tool [miloR](https://www.nature.com/articles/s41587-021-01033-z) allows you to look more deeply into subtle shifts in smaller neighborhoods of cells by utilizng differential abundance testing on the k-nearest neighbor graph.

<p align="center">
<img src="../img/milo_schematic.png" width="630">
</p>

# Creating Milo object

## Create new script

Next, open a new Rscript file, and start with some comments to indicate what this file is going to contain:

```r
# Single-cell RNA-seq analysis - differential abundance analysis with MiloR
```

Save the Rscript as `miloR_analysis_scrnaseq.R`.


## Load libraries

As usual, let us load the libraries needed at the beginning of our script.

```r
library(Seurat)
library(tidyverse)
library(SingleCellExperiment)
library(dplyr)
library(miloR)
```

## Select cell subsets

For continuity, let us take a look at the VSM cells and look at the differences between the `TN` and `cold7` conditions.

```r
seurat_vsm <- subset(seurat, subset = (celltype == "VSM"))
seurat_vsm <- subset(seurat_vsm, subset = (condition %in% c("TN", "cold7")))
```

MiloR generates the neighborhoods based upon the UMAP coordinates supplied, so we will re-run the necessary steps from our seurat pipeline on this new subset.

```r
# Re-calculate PCA and UMAP
seurat_vsm <- ScaleData(seurat_vsm)
seurat_vsm <- RunPCA(seurat_vsm, verbose = FALSE)
seurat_vsm <- RunUMAP(seurat_vsm, dims = 1:50)
seurat_vsm <- FindNeighbors(seurat_vsm, dims = 1:50)

Idents(seurat_vsm) <- "condition"
DimPlot(seurat_vsm)
```

<p align="center">
<img src="../img/miloR_umap_vsm.png" width="630">
</p>


## Creating single cell experiment

Seurat is not the only format with which we can load in single-cell data. There is another data structure known as `SingelCellExperiment` which MiloR makes use of. 

```r
sce_vsm <- as.SingleCellExperiment(seurat_vsm)
```

These objects have the following structure:

<p align="center">
<img src="../img/sce_description.png" width="630">
</p>

Image credit: (Amezquita, R.A., Lun, A.T.L., Becht, E. et al.)[https://doi-org.ezp-prod1.hul.harvard.edu/10.1038/s41592-019-0654-x]

We can use the functions from the SingleCellExperiment package to extract the different components. Let’s explore the counts and metadata for the experimental data.

```r
# Explore the raw counts for the dataset

## Check the assays present
assays(sce_vsm)

## Explore the raw counts for the dataset
dim(counts(sce_vsm))

counts(sce_vsm)[1:6, 1:6]
```

```
6 x 6 sparse Matrix of class "dgCMatrix"
       AAACCCACAGCTATTG_1 AAACGAAAGGGCGAAG_1 AAACGAACATTCGATG_1 AAACGAAGTAGCTCGC_1 AAACGAAGTTGGCTAT_1 AAACGAATCTGATTCT_1
Xkr4                    .                  .                  .                  .                  .                  .
Gm1992                  .                  .                  .                  .                  .                  .
Rp1                     .                  .                  .                  .                  .                  .
Sox17                   .                  .                  .                  .                  .                  .
Mrpl15                  2                  2                  .                  1                  1                  .
Lypla1                  .                  .                  .                  .                  .                  .
```

We see the raw counts data is a cell by gene sparse matrix with the same genes (rows) and columns (cells) as in our seurat object.

Next, we can get an idea of how to access the metadata in our object by using the `colData()` function:

```r
## Explore the cellular metadata for the dataset
dim(colData(sce_vsm))

head(colData(sce_vsm))
```

```
DataFrame with 6 rows and 17 columns
                    orig.ident nCount_RNA nFeature_RNA      sample log10GenesPerUMI mitoRatio   condition    S.Score   G2M.Score       Phase
                   <character>  <numeric>    <integer> <character>        <numeric> <numeric> <character>  <numeric>   <numeric> <character>
AAACCCACAGCTATTG_1   GSE160585      21744         4226    Sample_1         0.835980 0.0775846          TN -0.0357201 -0.04712786          G1
AAACGAAAGGGCGAAG_1   GSE160585       9727         3016    Sample_1         0.872506 0.0352590          TN -0.0435373 -0.04958206          G1
AAACGAACATTCGATG_1   GSE160585       7927         2319    Sample_1         0.863095 0.0586603          TN  0.0252971 -0.00971574           S
AAACGAAGTAGCTCGC_1   GSE160585      10858         2874    Sample_1         0.856963 0.0751520          TN -0.0367748 -0.02159758          G1
AAACGAAGTTGGCTAT_1   GSE160585       8131         2648    Sample_1         0.875394 0.0558357          TN -0.0354849 -0.00630089          G1
AAACGAATCTGATTCT_1   GSE160585       6899         2543    Sample_1         0.887089 0.0610233          TN -0.0366306 -0.02140319          G1
                   CC.Difference nCount_SCT nFeature_SCT integrated_snn_res.1.2    celltype seurat_clusters    ident
                       <numeric>  <numeric>    <integer>              <integer> <character>     <character> <factor>
AAACCCACAGCTATTG_1    0.01140775       6658         2106                      2         VSM               2       TN
AAACGAAAGGGCGAAG_1    0.00604473       7433         3006                      3         VSM               3       TN
AAACGAACATTCGATG_1    0.03501286       7175         2319                      3         VSM               3       TN
AAACGAAGTAGCTCGC_1   -0.01517724       7443         2823                      3         VSM               3       TN
AAACGAAGTTGGCTAT_1   -0.02918406       7264         2648                      3         VSM               3       TN
AAACGAATCTGATTCT_1   -0.01522739       6873         2543                      1         VSM               1       TN
```

## Creating Milo object

Now that we better understand how to use a SingleCellExperiment, we can convert it to a Milo object. While there are slight differences in this object, the basic idea of how to access metadata and counts information is consistent with a SingleCellExperiment.

```r
# Create miloR object
milo_vsm <- Milo(sce_vsm)
milo_vsm
```

```
class: Milo 
dim: 19771 5847 
metadata(0):
assays(2): counts logcounts
rownames(19771): Xkr4 Gm1992 ... CAAA01118383.1
  CAAA01147332.1
rowData names(0):
colnames(5847): AAACCCACAGCTATTG_1
  AAACGAAAGGGCGAAG_1 ... TTTGACTAGGCTTCCG_16
  TTTGTTGAGGGACAGG_16
colData names(17): orig.ident nCount_RNA ...
  seurat_clusters ident
reducedDimNames(0):
mainExpName: RNA
altExpNames(0):
nhoods dimensions(2): 1 1
nhoodCounts dimensions(2): 1 1
nhoodDistances dimension(1): 0
graph names(0):
nhoodIndex names(1): 0
nhoodExpression dimension(2): 1 1
nhoodReducedDim names(0):
nhoodGraph names(0):
nhoodAdjacency dimension(2): 1 1
```

# Milo workflow

Now that we have our dataset in the correct format, we can begin utilize the Milo workflow.


## Creating neighborhoods

Step one is to generate the k-nearest neighborhood graph with the `buildGraph()` function. The parameters include selected a `k` neighbors and `d` dimension values:

- `k`: An integer scalar that specifies the number of nearest-neighbours to consider for the graph building. Default is 10.
- `d`: The number of dimensions to use if the input is a matrix of cells. Deafult is 50.

Note that building the k-nearest neighbor graph may take some time.

```r
# Build the graph
traj_milo <- buildGraph(milo_vsm, k = 10, d = 30)
```

Afterwards we use the `makeNhoods()` function to define the neighborhoods themselves. To do so, the function randomly samples graph vertices, then refines them to collapse down the number of neighbourhoods to be tested.

- `prop`: A double scalar that defines what proportion of graph vertices to randomly sample. Must be 0 < `prop` < 1. Default is 0.1.
- `k`: An integer scalar - the same k used to construct the input graph. Default is 21.
- `d`: The number of dimensions to use if the input is a matrix of cells X reduced dimensions. Default is 30.

We additionally want to visualize how many cells belong to each neighborhood.

```r
traj_milo <- makeNhoods(traj_milo, prop = 0.1, k = 10, d=30, refined = TRUE)

plotNhoodSizeHist(traj_milo)
```

Now that we have identified which cells belong to which neighborhoods, we can quantify how many cells from each `sample` belong to each neighborhood.

```r
traj_milo <- countCells(traj_milo, 
                        meta.data = data.frame(colData(traj_milo)), 
                        sample="sample")

nhoodCounts(traj_milo) %>% head()
```

```
6 x 8 sparse Matrix of class "dgCMatrix"
  Sample_1 Sample_2 Sample_9 Sample_10 Sample_7 Sample_8 Sample_15 Sample_16
1       31        1        .         8        .        .         .         .
2        3        .       14        26        .        .         .         .
3        .        .        .        12        .        .         .         .
4       33        2        .         .        .        .         .         .
5        1        1        5        26        .        .         .         .
6        .        .        .         .        7        8        16         4
```

## Creating metadata

Next we need to create a metadata dataframe with all of the relevant pieces of information for the comparisons we want to run, including the sample names as well. In the case of this experiment, we need the columns `sample` and `condition`. 

```r
# Get columns of importance
traj_design <- data.frame(colData(traj_milo))[, c("sample", "condition")]

# Get unique values and set sample as the rownames
traj_design <- distinct(traj_design) %>% 
    remove_rownames() %>% 
    column_to_rownames("sample")

# Reorder rownames to match columns of nhoodCounts(milo)
traj_design <- traj_design[colnames(nhoodCounts(traj_milo)), , drop=FALSE]

traj_design
```

```
          condition
Sample_1         TN
Sample_2         TN
Sample_9         TN
Sample_10        TN
Sample_7      cold7
Sample_8      cold7
Sample_15     cold7
Sample_16     cold7
```

## Run differential abundance

Milo uses an adaptation of the Spatial FDR correction introduced by cydar, which accounts for the overlap between neighbourhoods. Specifically, each hypothesis test P-value is weighted by the reciprocal of the kth nearest neighbour distance. To use this statistic we first need to store the distances between nearest neighbors in the Milo object.

This calculates a Fold-change and corrected P-value for each neighbourhood, which indicates whether there is significant differential abundance between conditions.


```r
traj_milo <- calcNhoodDistance(traj_milo, d=30)
da_results <- testNhoods(traj_milo, 
                        design = ~condition, 
                        design.df = traj_design)

da_results %>% head()
```

```
      logFC   logCPM         F       PValue         FDR Nhood  SpatialFDR
1  5.383966 11.80962  3.997345 0.0456731668 0.150747056     1 0.147279465
2  5.439519 11.84264  4.145856 0.0418350534 0.150747056     2 0.147279465
3  3.514925 10.98646  1.607903 0.2048973695 0.252574320     3 0.249984729
4  5.353711 11.79050  3.950876 0.0469489658 0.150747056     4 0.147279465
5  5.063808 11.63182  3.576596 0.0587072653 0.150747056     5 0.147279465
6 -7.876050 13.52804 13.522319 0.0002403926 0.006379396     6 0.009554048
```

## Visualization

Lastly, we can visualize the differential abundance results using the UMAP coordinates that we initially supplied. 

Each dot represents a different neighborhood, with the size of the dot representing how many cells belong to that neighborhood. The color of the circle represent the log-fold change for that neighborhood, from the `da_results` object.

```r
traj_milo <- buildNhoodGraph(traj_milo)
plotNhoodGraphDA(traj_milo, da_results, alpha=0.05)
```

<p align="center">
<img src="../img/milo_umap.png" height="500">
</p>


---

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

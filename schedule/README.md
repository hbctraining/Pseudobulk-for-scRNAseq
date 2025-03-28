
# Workshop Schedule

## Pre-reading
* [Introduction to scRNA-seq](https://hbctraining.github.io/Intro-to-scRNAseq/lessons/01_intro_to_scRNA-seq.html)
* [scRNA-seq: From sequence reads to count matrix](https://hbctraining.github.io/Intro-to-scRNAseq/lessons/02_SC_generation_of_count_matrix.html)
* [scRNA-seq: From counts to clusters](../lessons/00_counts_to_clusters_overview.md)
* [Download this project](https://www.dropbox.com/scl/fi/ucxyecmviu4o882envynh/DGE_scRNA.zip?rlkey=ahb22wapzfkcbyioajlnq4m9k&st=sbfzxj8k&dl=1)

## Day 1

| Time |  Topic  | Instructor |
|:-----------:|:----------:|:--------:|
| 09:30 - 09:45 | [Workshop introduction](../lectures/workshop_intro_slides.pdf) |  Noor |
| 09:45 - 10:15 | Pre-reading review and Q&A | All |
| 10:15 - 10:25 | Break | |
| 10:25 - 11:00 | [Project setup and data exploration ](../lessons/01_setup_intro_dataset.md) | Noor |
| 11:00 - 11:50 | [Differential expression analysis using `FindMarkers()`](../lessons/02_DEanalysis_using_FindMarkers.md) | Meeta |
| 11:50 - 12:00 | Overview of self-learning materials and homework submission | Meeta |


### Before the next class:

I. Please **study the contents** and **work through all the code** within the following lessons:
 
   1. [Aggregating counts by celltype using pseudobulk approach](../lessons/03_pseudobulk_DESeq2.md)
      <details>
       <summary><i>Click here for a preview of this lesson</i></summary>
       <br> Forming pseudobulk samples is important to perform accurate differential expression analysis. Treating each cell as an independent replicate leads to underestimation of the variance and misleadingly small p-values. Working on the level of pseudobulk ensures reliable statistical tests because the samples correspond to the actual units of replication.  <br><br>In this lesson you will:<br>
             - Aggregate counts for a given celltype<br>
             - Demonstrate an efficent way to aggregate counts for multiple celltypes<br>
             - Use the aggregated counts to create a DESeq2 object for downstream analysis<br>
<br>
        </details>

  2. [DE analysis of pseudobulk data using DESeq2](../lessons/04_pseudobulk_DE_analysis.md)
      <details>
       <summary><i>Click here for a preview of this lesson</i></summary>
       <br> The next step is to take the DESeq2 object and run through the analysis workflow to identify differentially expressed genes. <br><br>In this lesson you will:<br>
           - Perform sample level QC<br>   
           - Evaluate gene-wise dispersions to evalute model fit<br>
           - Extract results and understand the statistics generated<br><br>
        </details>       
         

II. **Submit your work**:
   * Each lesson above contains exercises; please go through each of them.
   * **Submit your answers** to the exercises using [this Google form](https://forms.gle/sZeXdaUwf4uwwMah7) on **the day *before* the next class**.
   


### Questions?
* ***If you get stuck due to an error*** while running code in the lesson, [email us](mailto:hbctraining@hsph.harvard.edu) 

***

## Day 2

| Time |  Topic  | Instructor |
|:-----------:|:----------:|:--------:|
| 09:30 - 10:00 | Self-learning lessons discussion | All |
| 10:00 - 10:40 | [Visualization of differentially expressed genes](../lessons/05_pseudobulk_DE_visualizations.md) | Meeta |
| 10:40 - 10:50|  Break |  |
| 10:50 - 12:00 | [Comparison of results from different DE approaches](../lessons/06_DE_comparisons.md) | Noor|

### Before the next class:

I. Please **study the contents** and **work through all the code** within the following lessons:
   1. [Functional Analysis](../lessons/07_functional_analysis_pseudobulk.md)
      <details>
       <summary><i>Click here for a preview of this lesson</i></summary>
         <br>Now that we have significant genes, let's gain some biological insight <br><br>In this lesson, we will:<br>
             - Discuss approaches for functional analysis<br>
             - Use clusterProfiler to run over-representation analsyis and visualize results<br>
             - Use clusterProfiler to run GSEA <br><br>
        </details>

II. **Submit your work**:
   * Each lesson above contains exercises; please go through each of them.
   * **Submit your answers** to the exercises using [this Google form](https://forms.gle/7bmFiRZc6MQX7YXm7) on **the day *before* the next class**.

### Questions?
* ***If you get stuck due to an error*** while running code in the lesson, [email us](mailto:hbctraining@hsph.harvard.edu) 


***


## Day 3

| Time |  Topic  | Instructor |
|:-----------:|:----------:|:--------:|
| 09:30 - 10:15 | Self-learning lessons discussion | All |
| 10:15 - 11:15| [Methods for Differental Abundance](../lessons/08_differential_abundance.md)  | Noor  |
| 11:15 - 11:20 | Break |
| 11:25 - 12:00| Discussion and Q&A | All |
| 11:45 - 12:00| [Wrap-up](../lectures/Workshop_wrapup.pdf) | Meeta |

***

## Answer Keys
* [Day 1 exercises](https://hbctraining.github.io/DGE_analysis_scRNAseq/homework/Day_1_answer_key.html) 
* [Day 2 exercises](https://hbctraining.github.io/DGE_analysis_scRNAseq/homework/Day_2_answer_key.html) 





## Resources

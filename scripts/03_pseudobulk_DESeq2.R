# Single-cell RNA-seq analysis - Pseudobulk DE analysis with DESeq2


# Load libraries
library(Seurat)
library(tidyverse)
library(DESeq2)
library(pheatmap)
library(EnhancedVolcano)
library(RColorBrewer)
library(cowplot)
library(dplyr)


seurat <- readRDS("data/BAT_GSE160585_final.rds")

# Create metadata
meta_columns <- c("sample", "condition")
meta <- seurat@meta.data %>%
  select(meta_columns) %>%
  unique() %>%
  remove_rownames()
meta


# Aggregate count matrix by sample/celltype
bulk <- AggregateExpression(
  seurat,
  return.seurat = T,
  assays = "RNA",
  group.by = c("celltype", "sample", "condition")
)

# each 'cell' is a sample-condition-celltype pseudobulk profile
tail(Cells(bulk))


# Number of cells by sample and celltype
n_cells <- seurat@meta.data %>% 
  dplyr::count(sample, celltype) %>% 
  rename("n"="n_cells")
n_cells$sample <- str_replace(n_cells$sample, "_", "-")

# Add number of cells to bulk metadata
meta_bulk <- left_join(bulk@meta.data, n_cells)
rownames(meta_bulk) <- meta_bulk$orig.ident
bulk@meta.data <- meta_bulk

# Turn condition into a factor
bulk$condition <- factor(bulk$condition, levels=c("TN", "RT", "cold2", "cold7"))

bulk@meta.data %>% head()


bulk[["RNA"]]$counts[1:5, 1:5]


bulk_vsm <- subset(bulk, subset= (celltype == "VSM") & (condition %in% c("TN", "cold7")))


ggplot(bulk_vsm@meta.data, aes(x=sample, y=n_cells, fill=condition)) +
  geom_bar(stat="identity", color="black") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(x="Sample name", y="Number of cells") +
  geom_text(aes(label=n_cells), vjust=-0.5)


# Get count matrix
cluster_counts <- FetchData(bulk_vsm, layer="counts", vars=rownames(bulk_vsm))

# Create DESeq2 object
# transpose it to get genes as rows
dds <- DESeqDataSetFromMatrix(t(cluster_counts),
                              colData = bulk_vsm@meta.data,
                              design = ~ condition)

dds


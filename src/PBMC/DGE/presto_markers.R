#################################################
################### Set up ######################
#################################################

## load libraries

library(parallel)
library(gridExtra)
library(knitr)
library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)
library(ggExtra)
library(stats)
library(broom)
# library(RColorBrewer)
library(viridis)
library(scales)
library(ggrepel)
library(circlize)
library(RColorBrewer)
library(gplots)
library(cowplot)
library(gtable)
library(data.table)
library(Seurat)
library(hdf5r)
library(ggsci)
library(harmony)
library(future.apply)
library(speedglm)
library(rslurm)
library(lmtest)
library(presto)
library(here)





#########################################################
################### paths and stems #####################
#########################################################

# Get wd where rds are stored
setwd(here())


#########################################################
################### get seurat object ###################
#########################################################

results_dir <- "results/PBMC/seurat_objects/"
dir.create(results_dir, showWarnings = FALSE)

merged.seurat <- readRDS(file.path(results_dir,"curated_types_master_seurat_harmony_iterative_batch_emulsion.rds"))



# DGE on cluster vs all per main lineage
#--------------------------------------
# Select a cluster assignment iteration
cluster_assignments <- merged.seurat@meta.data$curated_types

# run subclustering on each cluster
presto_results <- lapply(unique(cluster_assignments), function(cluster) {
  
  # subset to given cluster
  tmp_seurat <- merged.seurat[, cluster_assignments == cluster]
  
  # run low dim embed functions
  
  general_results <- wilcoxauc(tmp_seurat, "curated_clusters", assay = 'data')
  
  
})

saveRDS(presto_results, file.path("results","DGE","presto_markers_curated_clusters.rds"))





#----------------------------
# DGE paired for each cluster


cluster_assignments <- merged.seurat@meta.data$curated_clusters

combinations <- combn(unique(cluster_assignments),
                      2,
                      simplify = F)

paired_results <- lapply(combinations, function(comb) {
  wilcoxauc(merged.seurat, "curated_clusters", groups_use = comb, assay = 'data') %>%
    mutate(groupA = comb[1]) %>%
    mutate(groupB = comb[2])
})

saveRDS(paired_results, file.path("results","DGE","presto_paired_markers_curated_clusters.rds"))





## For R3.5!


#################################################
################### Set up ######################
#################################################

## Take the input argument
# args <- commandArgs(TRUE)
# sample_name <- args[1]
# exp_matrix <- args[2]
# demux_best <- args[3]
# out_path <- args[4]

## load libraries

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
library(presto)



## ggplot theme 

theme_mrl <- function(x = 1) {
  theme_minimal() +
    theme(
      axis.line = element_line(),
      axis.ticks.x = element_line(),
      axis.ticks.y = element_line(),
      axis.text.x = element_text(size = 12*x,face = "bold", angle = 45, vjust = 0.6),
      axis.text.y = element_text(size = 12*x,face = "bold"),
      axis.title.x = element_text(size = 12*x,face = "bold"),
      axis.title.y = element_text(size = 12*x,face = "bold"),
      strip.background = element_rect(fill="gray20", colour="gray20", linetype="solid"),
      strip.text = element_text(size=14*x, colour="white", face="bold"),
      legend.title = element_text(size=14*x, face = "bold"),
      legend.text = element_text(size=12*x, color="gray20", face="bold"),
      legend.background = element_rect(fill = "transparent", colour = "transparent"),
      plot.title =  element_text(hjust=0.5, vjust=2, face="bold"),
      plot.subtitle = element_text(hjust=0.5, vjust=3, face="italic"),
      plot.caption = element_text(hjust = 0, face = "italic")
    )
}






#########################################################
################# Input/output info #####################
#########################################################

setwd(here::here())

#--------------------
# directory structure
results_dir <- "results/CSF_analysis"
dir.create(results_dir, showWarnings = FALSE)

objects_dir <- file.path(results_dir,"seurat_objects")
dir.create(objects_dir, showWarnings = FALSE)

clust_dir <- file.path(results_dir,"clustering_typing")
dir.create(clust_dir, showWarnings = FALSE)

#########################################################
###################### import data ######################
#########################################################


# get seurat object
merged.seurat <- readRDS(file.path(objects_dir,"merged.CSF_integrated_harmony_corrected_cluster_hcluster.rds"))

meta <- merged.seurat@meta.data %>% 
  tibble::rownames_to_column("index") %>%
  mutate(hcluster = replace_na(as.numeric(hcluster), "")) %>%
  mutate(cluster_label_h = paste(cluster_label, hcluster)) %>%
  separate(cluster_label_h, c("main"), sep = "_", extra = "drop", remove = F) %>%
  mutate(curated_clusters = unclass(forcats::fct_infreq(cluster_label_h))) %>%
  tibble::column_to_rownames("index")

merged.seurat@meta.data <- meta

############################
####### find markers ####### 
############################


#----------------------------
# general DGE for each cluster


presto_results <- wilcoxauc(merged.seurat, "curated_clusters", assay = 'data')

saveRDS(presto_results, file.path(clust_dir,"presto_markers_curated_clusters.rds"))



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

saveRDS(paired_results, file.path(clust_dir,"presto_paired_markers_curated_clusters.rds"))

paired_results <- do.call(rbind, paired_results) %>%
  dplyr::filter(padj < 0.01 & logFC > 0 & auc > 0.6)

saveRDS(paired_results, file.path(clust_dir,"presto_paired_markers_curated_clusters_filtered.rds"))


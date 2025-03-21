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


#--------------------
# input path
input_seurat <- file.path(objects_dir, "merged.CSF_integrated_harmony_corrected_cluster_hcluster.rds")

#--------------------
# output path

#output name stem
out_stem <- "CSF_extra"

# table output
output_table <- file.path(clust_dir,paste0("UMAP_harmony_",out_stem,".rds"))


#--------------------
# analysis options: 

# general
dims <- 1:30

#UMAP
umap_k <- c(20, 30, 50, 80)
umap_Md <- c(0.01,0.3,0.5,1,1.5, 2)
umap_spread <- c(0.1,0.5,1,3,6,10)

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
################### get seurat object ###################
#########################################################

merged.seurat <- readRDS(input_seurat)

cat("\n dim of seurat meta data \n")
dim(merged.seurat@meta.data)

results <- data.frame(cell.id = colnames(merged.seurat))


#########################################################
############## tSNE and UMAP reductions #################
#########################################################

cat("running tSNE")
merged.seurat <- RunTSNE(merged.seurat, reduction = "harmony", dims = dims)

results <- merged.seurat@reductions$tsne@cell.embeddings %>% 
  as.data.frame() %>%
  tibble::rownames_to_column("cell.id") %>%
  right_join(results)


cat("running UMAP iterations")

for(n in umap_k) {
  for(min.dist in umap_Md) {
    for(spread in umap_spread) {
      try(if(T){ 
        it <- paste0("UMAPn",n,"xMd",min.dist,"xS", spread)
        message(paste0("running ", it))
        merged.seurat <- RunUMAP(merged.seurat, 
                                 dims = dims, 
                                 reduction = "harmony", 
                                 n.neighbors = n,
                                 min.dist = min.dist, 
                                 spread = spread,
                                 reduction.key = it)
        
        results <- merged.seurat@reductions$umap@cell.embeddings %>% 
          as.data.frame() %>%
          tibble::rownames_to_column("cell.id") %>%
          left_join(results, by = "cell.id")
        cat("\n###########\nncol", ncol(results), "\n###########")
      })    		
    }
  }                        
}





#########################################################
#################### Save outputs #######################
#########################################################

saveRDS(results, output_table)


#qsave(merged.seurat, file.path(objects_dir,"master_seurat_harmony_low_res_clustered.qs"))

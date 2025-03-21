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
library(SingleR)
library(celldex)
library(BiocParallel)
library(here)

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
################### paths and stems #####################
#########################################################

# Get wd where rds are stored
setwd(here())



objects_dir <- "results/PBMC/seurat_objects/"
dir.create(objects_dir, showWarnings = FALSE)


results_dir <- "results/PBMC/clustering"
dir.create(results_dir, showWarnings = FALSE)

#########################################################
###################### import data ######################
#########################################################


# get seurat object
merged.seurat <- readRDS(file.path(objects_dir,"curated_types_master_seurat_harmony_iterative_batch_emulsion.rds"))



meta <- merged.seurat@meta.data %>%
    tibble::rownames_to_column("index") %>%
    mutate(curated_clusters = unclass(as.factor(meta_subC))) %>%
    tibble::column_to_rownames("index")
    
    
merged.seurat@meta.data <- meta



# Import singleR reference sets (needs internet connection, compute nodes have it on Ruddle)

monaco <- celldex::MonacoImmuneData()

dice <- celldex::DatabaseImmuneCellExpressionData()

#########################################################
################### Use SingleR refs ####################
#########################################################



# run singleR on the 2 ref sets
res.m <- SingleR(
            test = as.SingleCellExperiment(merged.seurat),
            ref = monaco,
            labels = monaco$label.fine)
            
res.d <- SingleR(
            test = as.SingleCellExperiment(merged.seurat),
            ref = dice,
            labels = dice$label.fine)

# convert to matrix and dataframes               
res.m <- list(scores = as.matrix(res.m$scores),
                calls = as.data.frame(res.m[,colnames(res.m) != "scores"]))
                
res.d <- list(scores = as.matrix(res.d$scores),
                calls = as.data.frame(res.d[,colnames(res.d) != "scores"]))


results <- list(res.m, res.d)


saveRDS(results, 
        file.path(results_dir,"cell_typing_sc_curated_types_master_seurat_harmony_iterative_batch_emulsion.rds"))

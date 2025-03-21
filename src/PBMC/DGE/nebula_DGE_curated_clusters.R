#################################################
################### Set up ######################
#################################################

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
library(scran)

#library(hdf5r)

library(ggsci)

library(harmony)
library(future.apply)

library(speedglm)
library(tradeSeq)
library(BiocParallel)
library(nebula)
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
################### nebula function #####################
#########################################################

nebula_wrapper <- function(count, preds, sid, formula, offset = c("lib_size", "scran", "none")) {

     # Decide offset method
  
      if(offset == "lib_size") {
        offset <- Matrix::colSums(count)
        offset <- offset / mean(offset)
      }
    
      if(offset == "scran") {
        offset <- scran::calculateSumFactors(count)
      }
      
      if(offset == "none") {
        offset <- NULL
      }


    # pred data for nebula formatted inputs
    
    sid <- pull(preds, sid)
    
    df <- model.matrix(formula, data=preds)
    
    # check grouping
    data_g = group_cell(count=count,id=sid,pred=df, offset = offset)
    
    # run nebula
    if(is.null(data_g)) {
        res <- nebula(count=count,id=sid,pred=df)
        
    } else {
        res <- nebula(data_g$count,data_g$id,pred=data_g$pred, offset = data_g$offset)
    }
    
    # Add avg gene expression to results
    avg_expr <- Matrix::rowSums(count) #sum count across cells
    avg_expr <- avg_expr / ncol(count) #divide by total number of cells
    avg_expr <- data.frame(avg_expr = avg_expr) %>% tibble::rownames_to_column("gene")
    
    
    res[['summary']] <- res[['summary']] %>% left_join(avg_expr)

    return(res)
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

MELD_LLH <- read_csv("results/PBMC/phate_meld_output/input_iterative_harmony_curated_clusters_MELD_LLH.csv")


# import updated metadata with module scores to filter potential contaminant cells
updated_meta <-readRDS(file.path(results_dir,"master_seurat_curated_types_table_harmony_iterative_batch_emulsion.rds"))

meta <- merged.seurat@meta.data %>%
    tibble::rownames_to_column("index") %>%
    left_join(MELD_LLH) %>%
    left_join(updated_meta) %>%
    mutate(curated_clusters = unclass(as.factor(meta_subC)))
    
cat("\n no metadata duplication:", ncol(merged.seurat) == nrow(meta))


#########################################################
################ Optional: remove cells #################
#########################################################

# Cluster 33 contaminants: moved to separate script
# meta <- meta %>%
#    mutate(cluster33_conta = curated_clusters == 33 & (naiveCD8 > -0.07 | IG > 0.005))

    
merged.seurat@meta.data <- meta %>%
    tibble::column_to_rownames("index")
    
# remove contaminants from object

# cat("\nnumber of cells initially", ncol(merged.seurat))

# merged.seurat <- merged.seurat[,!pull(meta,cluster33_conta)]

# cat("\nnumber of cells after removal", ncol(merged.seurat))

#########################################################
############## Choose clusters for loop #################
#########################################################


# Select a cluster assignment iteration
cluster_assignments <- merged.seurat@meta.data$curated_clusters

# print complete breakdown of all clusters as sanity check
meta %>% group_by(curated_clusters, meta_subC) %>% summarize(count = n()) %>% as.data.frame()


clusters_to_test <- meta %>% pull(curated_clusters) %>% unique()

cat("clusters_to_test:",clusters_to_test, sep = "\n")


#########################################################
##################### Run nebula ########################
#########################################################

# run nebula on each cluster
results <- lapply(clusters_to_test, function(cluster) {


    # subset to given cluster
    tmp_seurat <- merged.seurat[, cluster_assignments == cluster]
    
    cat("subset seurat for cluster", cluster, "to", ncol(tmp_seurat), "cells", "\n")
    
    # run low dim embed functions
    
    counts <- tmp_seurat@assays$RNA@counts
    
    preds <- tmp_seurat@meta.data
    
    nebula_wrapper(counts,preds, "sample", ~treatment+donor+nCount_RNA, offset = "none")
    
    
})




#########################################################
################# Save results to file ##################
#########################################################

names(results) <- paste0("curated_clusters.",clusters_to_test)

saveRDS(results, file.path("results",
                           "PBMC",
                           "DGE",
                           "nebula_results_allClusters_labelPred_countCov.rds"))



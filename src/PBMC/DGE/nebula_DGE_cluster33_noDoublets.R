# Initial analysis of CD16 monocytes (cluster 33) 
# showed some artefactual DGE results. This is likely due
# to some contaminants cells given the DEGs are IGs and typical 
# lymphoid lineage genes. Hence we do some filtering of
# contaminants (doublets) before rerunning NEBULA

#################################################
################### Set up ######################
#################################################

# whether to generate the contaminant filtering plots
plot=TRUE

## load libraries

library(gridExtra)
library(knitr)
library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)
library(magrittr)
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
        offset <- scran::computeSumFactors(count)
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
updated_meta <-readRDS(file.path(results_dir,"master_seurat_genesets_modulescores.rds"))

meta <- merged.seurat@meta.data %>%
    tibble::rownames_to_column("index") %>%
    left_join(MELD_LLH) %>%
    left_join(updated_meta) %>%
    mutate(curated_clusters = unclass(as.factor(meta_subC))) 
    
cat("\n no metadata duplication:", ncol(merged.seurat) == nrow(meta))




#########################################################
################ Optional: remove cells #################
#########################################################


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Prepare data to remove doublets from cluster 33 based on mixture of cell types

conta_genes <-
  c(
    grep(
      "^IGK|^IGL|^IGH|^TRBC|^TRBV|^TRBJ|^TRAC|^TRAV|^TRAJ",
      rownames(merged.seurat@assays$RNA@data),
      value = T
    ),
    c(
      "CD3E",
      "CD247",
      "GZMM",
      "MS4A1",
      "CD2",
      "CD8B",
      "ITK",
      "LCK",
      "GNLY"
    )
  )


# only keep cluster 33 cells
meta <- meta %>%
    mutate(conta_gene = Matrix::colSums(merged.seurat@assays$RNA@data[conta_genes,]) > 0) %>%
    dplyr::filter(curated_clusters == 33)

# import singleR results

mon_sc <- readRDS("results/PBMC/clustering/cell_typing_sc_curated_types_master_seurat_harmony_iterative_batch_emulsion.rds")[[1]]

# add barcodes to rownames in the score matrix
rownames(mon_sc[[1]]) <- rownames(mon_sc[[2]])

# only keep scores for cluster 33
cluster33_type_scores <- mon_sc[[1]][pull(meta, index),]

# add scores to the metadata
meta <- meta %>%
  left_join(cluster33_type_scores %>% as.data.frame() %>% tibble::rownames_to_column("index"))


# standardize for library size
cluster33_type_scores <- cluster33_type_scores / log10(meta$nCount_RNA)


# kmeans clustering on singleR scores

meta <- cluster33_type_scores %>%
  kmeans(2) %>%
  extract2("cluster") %>%
  data.frame(km = .) %>%
  tibble::rownames_to_column("index") %>%
  right_join(meta)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot the results

if(plot) {

p0 <- meta %>% 
  ggplot(aes(
    `Intermediate monocytes`,
    `Th1 cells`,
    color = as.factor(km)
  )) +
  geom_point() +
  scale_color_aaas(alpha = 0.5) +
  scale_x_log10() +
  theme_mrl(0.5) +
  theme(legend.position = "none")

  
p1 <- meta %>% 
  ggplot(aes(
    nCount_RNA,
    `Vd2 gd T cells` / log10(nCount_RNA),
    color = as.factor(km)
  )) +
  geom_point(alpha = 0.5) +
  scale_color_aaas() +
  scale_x_log10()+
  theme_mrl(0.5) +
  theme(legend.position = "none")


p2 <- meta %>% 
  ggplot(aes(
    nCount_RNA,
    `Th1 cells` / log10(nCount_RNA),
    color = as.factor(km)
  )) +
  geom_point() +
  theme_mrl(0.5) +
  scale_color_aaas(alpha = 0.5) +
  scale_x_log10()


p3 <- meta %>% 
  group_by(km) %>%
  mutate(km = paste0(km, " (", n(),")")) %>%
  ungroup() %>%
  ggplot(aes(as.factor(km),
    nCount_RNA,
    fill = as.factor(km)
  )) +
  geom_violin(alpha = 0.5) +
  scale_fill_aaas() +
  scale_y_log10() +
  theme_mrl(0.5) +
  theme(legend.position = "none")

p4 <- meta %>% 
  group_by(km) %>%
  mutate(km = paste0(km, " (", n(),")")) %>%
  ungroup() %>%
  ggplot(aes(as.factor(km),
    nFeature_RNA,
    fill = as.factor(km)
  )) +
  geom_violin(alpha = 0.5) +
  scale_fill_aaas() +
  scale_y_log10() +
  theme_mrl(0.5) +
  theme(legend.position = "none")

p5 <- meta %>% 
  group_by(km) %>%
  mutate(km = paste0(km, " (", n(),")")) %>%
  ungroup() %>%
  ggplot(aes(donor,
    t6M,
    fill = as.factor(km),
  )) +
  geom_violin(alpha = 0.5) +
  scale_fill_aaas() +
  facet_grid(.~km) +
  theme_mrl(0.5) +
  theme(legend.position = "none")

p6 <-
  plot_grid(plot_grid(p0,
                      p1,
                      p2 + theme(
                        legend.position = c(0, 1),
                        legend.justification = c(0, 1)
                      ),
                      nrow = 1),
            plot_grid(plot_grid(p3,
                                p4),
                      p5),
            nrow = 2)




ggsave(file.path("results", "PBMC", "DGE","cluster33_doublets.png"), p6)

}

# Cluster 33 contaminants:
meta <- meta %>%
    group_by(km) %>% 
    mutate(km = mean(nCount_RNA)) %>%
    ungroup() %>%
    mutate(km = km == min(km)) %>% # keep cells from cluster with lower lib size
    dplyr::filter(!conta_gene)
    
# remove contaminants from object

cat("\nnumber of cells initially", ncol(merged.seurat))

merged.seurat <- merged.seurat[,pull(meta,index)]

cat("\nnumber of cells after removal", ncol(merged.seurat))

    
merged.seurat@meta.data <- meta %>%
    tibble::column_to_rownames("index")
    
#stop("plotting done, stopping for control")

#########################################################
############## Choose clusters for loop #################
#########################################################


# Select a cluster assignment iteration
cluster_assignments <- merged.seurat@meta.data$curated_clusters

# print complete breakdown of all clusters as sanity check
meta %>% group_by(curated_clusters, meta_subC) %>% summarize(count = n()) %>% as.data.frame()


clusters_to_test <- 33

cat("clusters_to_test:",clusters_to_test, sep = "\n")


#########################################################
##################### Run nebula ########################
#########################################################

# run nebula on each cluster
results <- lapply(clusters_to_test, function(cluster) {


    # subset to given cluster
    tmp_seurat <- merged.seurat[, cluster_assignments == cluster]
    
    cat("subset seurat for cluster", cluster, "to", ncol(tmp_seurat), "cells", "\n")
    
    
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
                           "nebula_myeloid_results_v3_cluster33_NoContaminants_labelPred.rds"))



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



objects_dir <- "./results/PBMC/seurat_objects/"
dir.create(objects_dir, showWarnings = FALSE)


results_dir <- "./results/PBMC/clustering"
dir.create(results_dir, showWarnings = FALSE)

#########################################################
###################### import data ######################
#########################################################


# get seurat object
merged.seurat <- readRDS(file.path(objects_dir,"subclustered_master_seurat_harmony_iterative_batch_emulsion.rds"))


meta <- merged.seurat@meta.data %>% 
        tibble::rownames_to_column("cell.id")  %>%
        mutate(meta_subC = ifelse(metaclusters == "0", SubC_res0.2,
                                    ifelse(metaclusters == "CTLnCo", SubC_res0.2,
                                            ifelse(metaclusters == "1", SubC_res0.2,
                                                    ifelse(metaclusters == "4", SubC_res0.2,
                                                        ifelse(metaclusters == "3", SubC_res0.2,
                                                            ifelse(metaclusters == "myeloid", SubC_res0.05,
                                                                ifelse(metaclusters == "Bcells", SubC_res0.1,
                                                                    "unassigned" # should not happen, red flag is seen in final data
                                                                )
                                                            )
                                                        )
                                                    )
                                                )
                                            )
                                        )
                                    ) %>%
        mutate(curated_types = ifelse(meta_subC == "CTLnCo_2" | metaclusters == "1", "memCD4",
                                        ifelse(meta_subC %in% c("CTLnCo_6","myeloid_6"), "Bcells",
                                                ifelse(metaclusters == "3" | meta_subC == "myeloid_1", "CTLnCo",
                                                        ifelse(metaclusters == "0", "naiveCD4",
                                                                ifelse(metaclusters == "4", "naiveCD8",
                                                                    metaclusters
                                                                )
                                                        )
                                                )
                                        )
                            )
                ) %>%
        mutate(curated_clusters = unclass(meta_subC))






cat("number of cells in each curated type and clusters \n")
meta %>% group_by(curated_types,meta_subC) %>% summarize(count = n())

merged.seurat@meta.data <- meta %>% tibble::column_to_rownames("cell.id")




#########################################################
############## tSNE and UMAP reductions #################
#########################################################

# low embedding function

low_dim_embeddings <- function(seurat_object) {

    cat("running tSNE \n")
    seurat_object <- RunTSNE(seurat_object, reduction = "harmony", dims = 1:30)
    
    reductions.meta <- seurat_object@reductions$tsne@cell.embeddings %>% 
    		as.data.frame() %>%
    		tibble::rownames_to_column("cell.id")
    
    
    cat("running UMAP iterations")
    
    for(n in c(50)) { #knn param for UMAP
        for(min.dist in c(0.5)) {
            for(spread in c(5)) {
                try(if(T){ 
                    it <- paste0("UMAPn",n,"xMd",min.dist,"xS", spread)
                    message(paste0("running ", it))
                    seurat_object <- RunUMAP(seurat_object, 
                                            dims = 1:30, 
                                            reduction = "harmony", 
                                            n.neighbors = n,
                                            min.dist = min.dist, 
                                            spread = spread,
                                            reduction.key = it)
                                               
                    reductions.meta <- seurat_object@reductions$umap@cell.embeddings %>% 
                    		as.data.frame() %>%
                    		tibble::rownames_to_column("cell.id") %>%
                    		left_join(reductions.meta, by = "cell.id")
                })    		
            }
        }                        
    }
    #prefix columns to avoid confusion or conflicts with original results on whole set
    reductions.meta <- reductions.meta %>% tibble::column_to_rownames("cell.id")
    new_cols <- paste0("types_",colnames(reductions.meta))
    colnames(reductions.meta) <- new_cols
    
    
    
    return(reductions.meta %>% tibble::rownames_to_column("cell.id"))

}





#############################################################
########## run embeddings calc on main cell types ###########
#############################################################


# Select a cluster assignment iteration
cluster_assignments <- merged.seurat@meta.data$curated_types

# run subclustering on each cluster
curated_types_reductions <- lapply(unique(cluster_assignments), function(cluster) {

    # subset to given cluster
    tmp_seurat <- merged.seurat[, cluster_assignments == cluster]
    
    # run low dim embed functions
    
    reductions <- low_dim_embeddings(tmp_seurat)
    
    
})


# merge all the dataframe corresponding to each cluster
curated_types_reductions <- do.call(rbind, curated_types_reductions)


# add the results back to the main object
meta <- merged.seurat@meta.data %>% 
    tibble::rownames_to_column("cell.id") %>%
    left_join(curated_types_reductions, by = "cell.id") 

merged.seurat@meta.data <- meta %>%
    tibble::column_to_rownames("cell.id")
    


#########################################################
#################### Save outputs #######################
#########################################################

saveRDS(meta, file.path(results_dir,"master_seurat_curated_types_table_harmony_iterative_batch_emulsion.rds"))


saveRDS(merged.seurat, file.path(objects_dir,"curated_types_master_seurat_harmony_iterative_batch_emulsion.rds"))



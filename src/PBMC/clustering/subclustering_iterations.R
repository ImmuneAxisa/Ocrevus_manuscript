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
################### get seurat object ###################
#########################################################

merged.seurat <- readRDS(file.path(objects_dir,"low_res_clustered_master_seurat_harmony_iterative_batch_emulsion.rds"))

# metaclusters for main cell types:

meta <- merged.seurat@meta.data %>% 
        tibble::rownames_to_column("cell.id")  %>%
        mutate(treatment = ifelse(grepl("6", treatment), "t6M","t0M")) %>%
        mutate(SNNk20_algo1_res0.2 = as.numeric(as.character(SNNk20_algo1_res0.2))) %>%
        mutate(metaclusters = ifelse(SNNk20_algo1_res0.2 %in% c(2,7,9,10,11,14,15),
                                        "myeloid",
                                        ifelse(SNNk20_algo1_res0.2 %in% c(5,13),
                                                "Bcells",
                                                ifelse(SNNk20_algo1_res0.2 %in% c(6,8,12),
                                                        "CTLnCo",
                                                        SNNk20_algo1_res0.2))))
                                                        

cat("number of cells in each metaclusters \n")
summary(as.factor(meta$metaclusters))

merged.seurat@meta.data <- meta %>% tibble::column_to_rownames("cell.id")

# for testing
#merged.seurat <- merged.seurat[,1:1000]

#########################################################
########### clustering and cluster markers ##############
#########################################################

print("running clustering")


# clustering function

run_clustering <- function(seurat_object){

    print("find neighbors")
    for(k in c(20)) {
    
            seurat_object <- FindNeighbors(seurat_object,
                                            dims = 1:30, 
                                            k.param = k,
                                            reduction = "harmony")
    
    
    print("find clusters")
        for(algo in c(1)) {
                        
                        
                        # Loop with different resolutions to define clusters
                        
                        for(res in c(0.01,0.05,0.1,0.2)) {
                        
                        
                            seurat_object <- FindClusters(seurat_object, 
                                                            resolution = res,
                                                            algorithm = algo)
                            clusters <- seurat_object[["seurat_clusters"]][,1]
                            seurat_object[[paste0("SNNk",k,"_algo",algo,"_res",res)]] <- clusters
                            }
                            
                        }
                        
    }


    return(seurat_object)
}

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
        for(min.dist in c(0.01,0.3,0.5)) {
            for(spread in c(5,10)) {
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
    new_cols <- paste0("SubC_",colnames(reductions.meta))
    colnames(reductions.meta) <- new_cols
    
    
    
    return(reductions.meta %>% tibble::rownames_to_column("cell.id"))

}


#############################################################
###### run subclustering on low resolution communities ######
#############################################################


# Select a cluster assignment iteration
cluster_assignments <- merged.seurat@meta.data$metaclusters

# run subclustering on each cluster
subclustering <- lapply(unique(cluster_assignments), function(cluster) {


    # subset to given cluster
    tmp_seurat <- merged.seurat[, cluster_assignments == cluster]
    
    # rerun SNN graph building
    tmp_seurat <- FindNeighbors(tmp_seurat,
                                    dims = 1:30, 
                                    k.param = 20,
                                    reduction = "harmony")

    # find clusters at different resolutions
    for(res in c(0.05,0.1,0.2)) {
                        
                        
        tmp_seurat <- FindClusters(tmp_seurat, 
                                        resolution = res,
                                        algorithm = 1)
        
        # extract clustering results                                
        subclusters <- tmp_seurat[["seurat_clusters"]][,1]
        
        # rename it and add back to metadata
        tmp_seurat[[paste0("SubC_res",res)]] <- paste0(cluster,"_",subclusters)
    }
    
    
    # run low dim embed functions
    
    reductions <- low_dim_embeddings(tmp_seurat)
    
    # only keep the columns with the results of the subclustering, 
    # merge with low dim embed results,
    # and return the df on subcluster calls
    tmp_meta <- tmp_seurat@meta.data %>%
        tibble::rownames_to_column("cell.id") %>%
        dplyr::select(matches("SubC_res|cell.id")) %>%
        left_join(reductions, by = "cell.id")
    
})


# merge all the dataframe corresponding to each cluster
subclustering <- do.call(rbind, subclustering)


# add the results back to the main object
meta <- merged.seurat@meta.data %>% 
    tibble::rownames_to_column("cell.id") %>%
    left_join(subclustering, by = "cell.id") 

merged.seurat@meta.data <- meta %>%
    tibble::column_to_rownames("cell.id")


#########################################################
#################### Save outputs #######################
#########################################################

saveRDS(meta, file.path(results_dir,"master_seurat_subclustering_table_harmony_iterative_batch_emulsion.rds"))


saveRDS(merged.seurat, file.path(objects_dir,"subclustered_master_seurat_harmony_iterative_batch_emulsion.rds"))

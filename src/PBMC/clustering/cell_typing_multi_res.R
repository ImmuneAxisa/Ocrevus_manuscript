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
library(SingleR)
library(celldex)
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
setwd(here)



objects_dir <- "results/PBMC/seurat_objects/"
dir.create(objects_dir, showWarnings = FALSE)


results_dir <- "results/PBMC/clustering"
dir.create(results_dir, showWarnings = FALSE)

#########################################################
###################### import data ######################
#########################################################


merged.seurat <- readRDS(file.path(objects_dir,"subclustered_master_seurat_harmony_iterative_batch_emulsion.rds"))

original_meta <- merged.seurat@meta.data %>% 
        tibble::rownames_to_column("cell.id")  %>%
        dplyr::select(-metaclusters)


# get seurat object
merged.seurat <- readRDS(file.path(objects_dir,"subclustered_master_seurat_harmony_iterative_batch_emulsion_v2.rds"))


meta <- merged.seurat@meta.data %>% 
        tibble::rownames_to_column("cell.id")  %>%
        dplyr::select(cell.id,metaclusters, SubC_res0.5, SubC_res0.01) %>%
        left_join(original_meta) %>%
        mutate(meta_subC = ifelse(metaclusters == "0", SubC_res0.2,
                                    ifelse(metaclusters == "CTLnCo", SubC_res0.5,
                                            ifelse(metaclusters == "1", SubC_res0.2,
                                                    ifelse(metaclusters == "4", SubC_res0.2,
                                                        ifelse(metaclusters == "3", SubC_res0.2, #not present in v2 so should not match anything
                                                            ifelse(metaclusters == "myeloid", SubC_res0.01,
                                                                ifelse(metaclusters == "Bcells", SubC_res0.1,
                                                                    "unassigned"
                                                                )
                                                            )
                                                        )
                                                    )
                                                )
                                            )
                                        )
                                    )

cat("summary of metadata table")
summary(meta)



cat("number of cells in each meta subclusters \n")
summary(as.factor(meta$meta_subC))

merged.seurat@meta.data <- meta %>% tibble::column_to_rownames("cell.id")


# specify cluster call you want to use for cell typing with singleR 
cluster_calls <- c("SubC_res0.05", "SubC_res0.1", "SubC_res0.2", "meta_subC")

# Import singleR reference sets (needs internet connection, compute nodes have it on Ruddle)

monaco <- celldex::MonacoImmuneData()

dice <- celldex::DatabaseImmuneCellExpressionData()

#########################################################
################### Use SingleR refs ####################
#########################################################


results <- lapply(cluster_calls, function(cluster_index) {
    
    
        # run singleR on the 2 ref sets
        res.m <- SingleR(
                    test = as.SingleCellExperiment(merged.seurat),
                    ref = monaco,
                    labels = monaco$label.fine,
                    clusters = merged.seurat@meta.data %>% pull(cluster_index))
                    
        res.d <- SingleR(
                    test = as.SingleCellExperiment(merged.seurat),
                    ref = dice,
                    labels = dice$label.fine,
                    clusters = merged.seurat@meta.data %>% pull(cluster_index))
    
        # convert to matrix and dataframes               
        res.m <- list(scores = as.matrix(res.m$scores),
                        calls = as.data.frame(res.m[,colnames(res.m) != "scores"]))
                        
        res.d <- list(scores = as.matrix(res.d$scores),
                        calls = as.data.frame(res.d[,colnames(res.d) != "scores"]))


    list(res.m, res.d)
    
    
    
})

names(results) <- cluster_calls


saveRDS(results, 
        file.path(results_dir,"subclustered_master_seurat_harmony_iterative_batch_emulsion_singleR_v2.rds"))






#################################################
################### Set up ######################
#################################################

## load libraries

library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(data.table)
library(here)

#########################################################
################### paths and stems #####################
#########################################################


# Get wd where rds are stored
setwd(here())

#########################################################
############ write pca and metadata table  ##############
#########################################################

#-----------
# import data
#-----------
merged.seurat <- readRDS("results/CSF_analysis/seurat_objects/merged.CSF_integrated_harmony_corrected_cluster_hcluster.rds")

#---------------------------
# merge pc loadings and meta
#---------------------------
master_harmony <- merged.seurat@reductions$harmony@cell.embeddings %>% 
		as.data.frame() %>%
		tibble::rownames_to_column("cell.id") %>%
		left_join(merged.seurat@meta.data %>% tibble::rownames_to_column("cell.id"), by = "cell.id") %>%
        dplyr::rename("index" = "cell.id")
        
cat("\nnumber of cells initially:", 
    nrow(master_harmony),
    "\n\n")

master_harmony <- master_harmony %>%
        mutate(hcluster = as.character(hcluster)) %>%
        mutate(hcluster = replace_na(hcluster, "")) #%>%
        #dplyr::filter(hcluster != "1")

cat("\nnumber of cells after removing T cells in Mo cluster:", 
    nrow(master_harmony),
    "\n\n")

#-----------------------------        
# define splitting categories
#-----------------------------

cats <- unique(master_harmony$LP)

combinations <- combn(cats, 2, simplify = F)

print(combinations)

for(combinations_index in combinations) {


    # do filterting and cleaning before export for phate meld vfc        
            
    harmony <- master_harmony %>%
            dplyr::filter(LP %in% combinations_index)
    
    
    # split and write to file
      
    harmony_PCA <- harmony[,1:31]
    
    harmony_meta <- harmony[,c(1,52:ncol(harmony))]
    
    basename <- paste(combinations_index, collapse = "vs")
    
    
    write_csv(harmony_PCA, paste0("results/CSF_analysis/phate_meld_input_csv/",basename,"_pca.csv"))
    
    write_csv(harmony_meta, paste0("results/CSF_analysis/phate_meld_input_csv/",basename,"_meta.csv"))
    
    
    # run the shell wrapper to submit phate and meld jobs
    
    
    system2(command = "/gpfs/ycga/project/hafler/pa326/Ocrevus/src/CSF_analysis/phate_meld/phate_meld_sbatcher_w_args_v2.sh", 
            args    = c(here(), #project dir
                        here("results/CSF_analysis"), # results dir for CSF,
                        basename,
                        "LP", #column to use in meta to run MELD with
                        "cluster_label", #column to use in meta to split data into for the VFC 
                        paste0(basename, "_all_CSF"), #suffix to add for the output
                        
                        ))


}




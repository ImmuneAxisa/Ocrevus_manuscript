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

merged.seurat <- readRDS("results/PBMC/seurat_objects/curated_types_master_seurat_harmony_iterative_batch_emulsion.rds")



# merge pc loadings and meta

harmony <- merged.seurat@reductions$harmony@cell.embeddings %>% 
		as.data.frame() %>%
		tibble::rownames_to_column("cell.id") %>%
		left_join(merged.seurat@meta.data %>% tibble::rownames_to_column("cell.id"), by = "cell.id") %>%
        dplyr::rename("index" = "cell.id") %>%
        mutate(treatment = ifelse(grepl("6", treatment), "t6M","t0M")) %>%
        mutate(curated_clusters = unclass(as.factor(meta_subC)))
  
cat("number of cells in each clusters for VFC")
summary(as.factor(harmony$curated_clusters))

# split and write to file
  
harmony_PCA <- harmony[,1:31]

harmony_meta <- harmony[,c(1,52:ncol(harmony))]


write_csv(harmony_PCA, "results/PBMC/phate_meld_input_csv/pbmc_pca.csv")

write_csv(harmony_meta, "results/PBMC/phate_meld_input_csv/pbmc_meta.csv")


# run the shell wrapper to submit phate and meld jobs


system2(command = here("src/wrapper_scripts/phate_meld_sbatcher_w_args.sh"), 
        args    = c(here(), #project dir
                    here("results/PBMC"), # results dir for PBMC
                    "pbmc", #iteration basenames
                    "treatment", #column to use in meta to run MELD with
                    "curated_clusters", #column to use in meta to split data into for the VFC 
                    "iterative_harmony_curated_clusters" #suffix to add for the output
                    ))



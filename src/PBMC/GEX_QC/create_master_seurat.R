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
library(here)

setwd(here())



#################################################
################ Import data ####################
#################################################

# Set path to cell ranger outputs  

path <- "/gpfs/ycga/scratch60/hafler/pa326/Ocrevus/PBMC/"

samples <- list.files(path)

## Import expression matrix as a list of seurat objects
master_seurat <- lapply(samples, function(x) {
                h5.path <- paste0(path,x,"/filtered_feature_bc_matrix.h5")
                Seurat_object <- CreateSeuratObject(Read10X_h5(h5.path))
                return(Seurat_object)
})

# name list elements with sample names
names(master_seurat) <- gsub("_HHT_cellranger","", samples)


#################################################
############### merge and norm ##################
#################################################

# merge all seurat objects and add a suffix to each cell barcode with the sample names it belongs to
merged.seurat <- merge(master_seurat[[1]], master_seurat[-1], add.cell.ids = names(master_seurat))


saveRDS(merged.seurat, "results/PBMC/seurat_objects/master_seurat.rds")

saveRDS(merged.seurat@meta.data, "results/PBMC/seurat_QC/master_meta.rds")



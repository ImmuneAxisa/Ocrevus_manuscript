
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
library(Seurat)
library(DESeq2)
library(forcats)
library(here)
library(tictoc)


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
setwd(here::here())

# path to save files to
results_dir <- "results/CSF_analysis/DGE/poisson_labelPred_rbCov/"
dir.create(results_dir, showWarnings = FALSE, recursive = T)



#########################################################
############## import data and metadata  ################
#########################################################

merged.seurat <- readRDS("results/CSF_analysis/seurat_objects/merged.CSF_integrated_harmony_corrected_cluster_hcluster.rds")

meta <- merged.seurat@meta.data

GEX <- merged.seurat@assays$RNA@counts

rm(merged.seurat)
gc()

#########################################################
############ Reorganize cluster metadata ################
#########################################################


meta <- meta %>% 
    tibble::rownames_to_column("index") %>%
    mutate(hcluster = as.character(hcluster)) %>%
    mutate(hcluster = replace_na(hcluster, "")) %>%
    mutate(cluster_label = paste0(cluster_label,hcluster)) %>% 
    mutate(cluster_label = gsub("Mø","Mo",cluster_label))
    
meta <- meta %>% 
    separate(cluster_label, "main", sep = "_", remove = F, extra = "drop")
    
#meta$percent.rb <- PercentageFeatureSet(merged.seurat, pattern = "^RP[SL]")

cat("\nnumber of cells initially:", 
    nrow(meta),
    "\n\n")
    


############################################################################ 
############################################################################
############## Prepare DGE analysis ########################################
############################################################################
############################################################################ 

deseq2_dge_wrapper <- function(meta, GEX_mat) {
  
  cat("processing", nrow(meta), "samples\n\n")
  tic(quiet = FALSE)
  
  # select cells from cluster 
  GEX_mat <- GEX_mat[,meta$index]
  # remove low expression genes
  GEX_mat <- GEX_mat[Matrix::rowMeans(GEX_mat) > 0.005,]
  
  cat("testing", nrow(GEX_mat), "features\n\n")
  # build DDS object
  dds_global <- meta %>%
    mutate(disease = Donor_Type,
           Trt = ifelse(Donor_Type == "Healthy", "Baseline", LP)) %>%
    group_by(disease) %>%
    mutate(patient = paste0("ID",unclass(as.factor(patient)))) %>%
    dplyr::select(index, disease, Trt, patient) %>%
    tibble::column_to_rownames("index") %>%
    DESeqDataSetFromMatrix(as.matrix(GEX_mat)[,rownames(.)], .,
                           ~1)
  
  #-------------------------------
  # build model matrix beforehand
  # to manually drop unused levels
  #-------------------------------
  mm <- model.matrix(~disease+disease:patient+disease:Trt, colData(dds_global))
  
  all.zero <- apply(mm, 2, function(x) all(x==0))
  idx <- which(all.zero)
  mm <- mm[,-idx]
  
  design(dds_global) <- mm
  
  colnames(mm)
  
  #--------------------------------------------
  # run full model to extract expression values
  #--------------------------------------------
  full_model <- DESeq(dds_global, fitType = "glmGamPoi",
                      full = mm)
  
  #-------------------------------
  # get pre-post 
  #-------------------------------
  # as recommended we use the LRT test to get significance 
  results_pre_post_global <- DESeq(dds_global, fitType = "glmGamPoi", test = "LRT", 
                                   full = mm, reduced = mm[,colnames(mm) != "diseaseMS:TrtFollow_up"]) %>%
    results(independentFiltering = TRUE) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene_name")
  
  # but we also at least want the size effect estimates
  results_pre_post_global <- full_model %>%
    lfcShrink(coef = "diseaseMS.TrtFollow_up", type = "ashr") %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene_name") %>%
    dplyr::select(gene_name, slfc = log2FoldChange) %>%
    right_join(results_pre_post_global)
  
  
  #-------------------------------
  # get disease genes
  #-------------------------------
  # as recommended we use the LRT test to get significance 
  results_disease_global <- DESeq(dds_global, fitType = "glmGamPoi", test = "LRT", 
                                  full = mm, reduced = mm[,colnames(mm) != "diseaseMS"]) %>%
    results(independentFiltering = TRUE) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene_name")
  
  # but we also at least want the size effect estimates
  results_disease_global <- full_model %>%
    lfcShrink(coef = "diseaseMS", type = "ashr") %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene_name") %>%
    dplyr::select(gene_name, slfc = log2FoldChange) %>%
    right_join(results_disease_global)
  
  
  rm(full_model)
  gc()
  toc()
  
  return(list(disease = results_disease_global,
              Trt = results_pre_post_global))
  
}


############

meta %>%
  group_by(cluster_label) %>%
  nest() %>%
  mutate(dge_results = purrr::map(data, 
                                  deseq2_dge_wrapper, 
                                  GEX_mat = GEX, 
                                  .progress = TRUE)) %>%
  qs::qsave("reviewing/deseq2_DGE_CSF.qs")


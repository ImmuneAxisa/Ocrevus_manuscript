# Prepare cleaned up object with relevant metadata from the analysis



library(tidyverse)
library(Seurat)
library(ggsci)
library(colorspace)


################################################################################
################################## PBMC data ###################################
################################################################################


source("src/wrapper_scripts/results_annotations.R")

best_UMAP <- best_UMAP %>%
  dplyr::select(cell.id, UMAP_1 = dim_1, UMAP_2 = dim_2)

meld <- read_csv("results/PBMC/phate_meld_output/input_iterative_harmony_curated_clusters_MELD_LLH.csv") %>%
  dplyr::select(cell.id = index, meld_t6M = t6M)

meta_table <- meta_table %>%
  dplyr::select(-matches("^coarse|^Sub|^SNN|^RNA|^meta|UMAP|tSNE"), -c(barcode, orig.ident,seurat_clusters)) %>%
  left_join(best_UMAP) %>%
  left_join(cluster_labels_table) %>%
  left_join(meld)

tibble::column_to_rownames(meta_table,"cell.id") %>%
  colnames() %>%
  paste0(":") %>%
  cat("PBMC meta data column names:\n", ., file = "meta_data_dictionary.txt", sep = "\n")


merged.seurat <- readRDS("results/PBMC/seurat_objects/master_filtered_seurat.rds")

merged.seurat@meta.data <- tibble::column_to_rownames(meta_table,"cell.id")[colnames(merged.seurat),]

saveRDS(merged.seurat, "PBMC_seurat.rds")



################################################################################
################################## CSF data ####################################
################################################################################

seurat <- readRDS("results/CSF_analysis/seurat_objects/merged.CSF_integrated_harmony_corrected_cluster_hcluster.rds")

new_names <- tribble(
  ~cluster_label_h, ~cluster_label_new,
  "B cell_1 ", "memory B cells",
  "B cell_2 ", "plasmablasts",
  "CD4_1 ", "memory CD4 1",
  "CD4_2 ", "memory CD4 2",
  "CD4_3 ", "memory CD4 3",
  "CD4_4 ", "memory CD4 4",
  "CD8_1 ", "memory CD8 1",
  "CD8_2 ", "memory CD8 2",
  "DC_1 ", "DC_1",
  "DC_2 ", "DC_2",
  "DC_3 ", "DC_3",
  "gd T ", "gd T",
  "Mø 1", "memory CD4 5",
  "Mø 2", "Mac_1",
  "Monocytes ", "Mac_2",
  "NK ", "NK",
  "pDC ", "pDC",
  "Treg ", "Treg"
)

umaps_cleaned <- readRDS("results/CSF_analysis/clustering_typing/UMAP_harmony_CSF_extra.rds") %>%
  dplyr::select(matches("cell\\.id|^UMAPn30")) %>%
  pivot_longer(-cell.id, names_to = "it", values_to = "value") %>%
  separate(it, c("it", "dim"), sep = "_") %>%
  mutate(dim = paste0("dim",dim)) %>%
  pivot_wider(names_from = dim, values_from = value) %>%
  dplyr::filter(!grepl("tSNE", it)) %>%
  group_by(it) %>%
  mutate(dim1 = scales::rescale(dim1), dim2 = scales::rescale(dim2)) %>%
  separate(it, c("it", "md", "s"), sep = "x") %>%
  mutate(md = parse_number(gsub("Md([0-9])", "\\1\\.", md)),
         s = parse_number(s)) %>%
  dplyr::filter(it == "UMAPn30" & md == 2 & s == 1) %>%
  dplyr::select(cell.id, UMAP_1 = dim1, UMAP_2 = dim2)

BvF <- read_csv("./results/CSF_analysis/phate_meld_output/BaselinevsFollow_up_BaselinevsFollow_up_all_CSF_MELD_LLH.csv") %>%
  dplyr::select(index, meld_Follow_up = Follow_up)
BvHC <- read_csv("results/CSF_analysis/phate_meld_output/BaselinevsHealthy_BaselinevsHealthy_all_CSF_MELD_LLH.csv") %>%
  dplyr::select(index, meld_MS_preTx = Baseline)

meta <- seurat@meta.data %>% 
  tibble::rownames_to_column("index") %>%
  mutate(hcluster = replace_na(as.character(hcluster), "")) %>%
  mutate(cluster_label_h = paste(cluster_label, hcluster)) %>%
  separate(cluster_label_h, c("main"), sep = "_", extra = "drop", remove = F) %>%
  mutate(curated_clusters = unclass(forcats::fct_infreq(cluster_label_h))) %>%
  left_join(new_names) %>%
  dplyr::select(-matches("^RNA|^celltype"), -c(hcluster, cluster_label_h, SingleR.cluster.labels, cluster_label, seurat_clusters)) %>%
  dplyr::rename(main_lineage = main, cluster_label = cluster_label_new) %>%
  left_join(umaps_cleaned, by = c("index" = "cell.id")) %>%
  left_join(BvF) %>%
  left_join(BvHC) %>%
  tibble::column_to_rownames("index") 

seurat@meta.data <- meta[colnames(seurat),]


saveRDS(seurat, "CSF_seurat.rds")

meta %>%
  colnames() %>%
  paste0(":") %>%
  cat("\nCSF meta data column names:\n", ., file = "meta_data_dictionary.txt", sep = "\n", append = TRUE)

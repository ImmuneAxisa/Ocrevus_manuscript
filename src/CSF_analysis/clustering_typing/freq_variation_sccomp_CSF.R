library(tidyverse)
library(ggsci)
library(colorspace)

setwd(here::here())

source("src/wrapper_scripts/results_annotations.R")

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
  "Mø 2", "CSF macrophages",
  "Monocytes ", "monocytes",
  "NK ", "NK",
  "pDC ", "pDC",
  "Treg ", "Treg"
)

cluster_colors <- alpha(pals::cols25(), 0.8)[1:nrow(new_names)]

names(cluster_colors) <- pull(new_names, cluster_label_new)

cluster_colors <- cluster_colors[order(names(cluster_colors))]



seurat <- readRDS("results/CSF_analysis/seurat_objects/merged.CSF_integrated_harmony_corrected_cluster_hcluster.rds")

meta <- merged.seurat@meta.data %>% 
  tibble::rownames_to_column("index") %>%
  mutate(hcluster = replace_na(as.character(hcluster), "")) %>%
  mutate(cluster_label_h = paste(cluster_label, hcluster)) %>%
  separate(cluster_label_h, c("main"), sep = "_", extra = "drop", remove = F) %>%
  mutate(curated_clusters = unclass(forcats::fct_infreq(cluster_label_h))) %>%
  left_join(new_names) %>%
  unite("ID", patient, LP, sep = "___", remove = FALSE)




library(sccomp)

sccomp_disease <- meta %>%
  dplyr::filter(LP != "Follow_up") %>%
  sccomp_glm(
    formula_composition = ~ LP,
    .sample = ID,
    .cell_group = cluster_label_new,
    bimodal_mean_variability_association = TRUE,
    cores = 1
  )

sccomp_trt <- meta %>%
  dplyr::filter(LP != "Healthy") %>%
  sccomp_glm(
    formula_composition = ~ LP+patient,
    .sample = ID,
    .cell_group = cluster_label_new,
    bimodal_mean_variability_association = TRUE,
    cores = 1
  )

saveRDS(list(trt = sccomp_trt, disease = sccomp_disease), 
        "reviewing/sccomp_paired_CSF.rds")

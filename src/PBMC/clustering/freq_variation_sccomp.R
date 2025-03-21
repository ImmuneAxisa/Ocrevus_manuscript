library(tidyverse)
library(ggsci)
library(colorspace)

setwd(here::here())

source("src/wrapper_scripts/results_annotations.R")

cluster_groups <- meta_table %>%
  dplyr::select(cell.id, coarse_clusters, curated_clusters, donor, treatment, sample) %>%
  left_join(cluster_labels_table)



library(sccomp)

sccomp_paired <- sccomp_glm(cluster_groups,
                            formula_composition = ~treatment+donor,
                            .sample = sample,
                            .cell_group = fine_clusters,
                            bimodal_mean_variability_association = TRUE,
                            cores = 1)

saveRDS(sccomp_paired, "reviewing/sccomp_paired.rds")

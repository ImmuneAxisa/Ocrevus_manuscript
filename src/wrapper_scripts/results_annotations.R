# annotation helpers



meta_table <- readRDS("results/PBMC/clustering/master_seurat_curated_types_table_harmony_iterative_batch_emulsion.rds") %>%
  mutate(treatment = ifelse(grepl("6", treatment), "t6M","t0M"),
         curated_clusters = unclass(as.factor(curated_clusters)))



meta_table <- readRDS("results/PBMC/clustering/master_seurat_coarse_types_table_harmony_iterative_batch_emulsion.rds") %>%
  dplyr::select(matches("coarse|cell\\.id")) %>%
  right_join(meta_table)



best_UMAP <- readRDS("results/PBMC/clustering/best_UMAP.rds")


cluster_labels_table <- tibble(
  "granularity" = c("coarse_clusters", "fine_clusters"),
  "1" = c("naiveT", "naiveCD4 1"),
  "2" = c("naiveT", "naiveCD4 2"),
  "3" = c("naiveT", "naiveCD4 3"),
  "4" = c("naiveT", "naiveCD4 4"),
  "5" = c("memCD4", "memCD4 1"),
  "6" = c("memCD4", "memCD4 2"),
  "7" = c("memCD4", "Treg"),
  "8" = c("memCD4", "memCD4 3"),
  "9" = c("memCD4", "memCD4 4"),
  "10" = c("memCD4", "memCD4 5"),
  "11" = c("Cytotoxic", "memCD8 1"),
  "12" = c("Cytotoxic", "memCD8 2"),
  "13" = c("Cytotoxic", "MAIT"),
  "14" = c("Cytotoxic", "gdT"),
  "15" = c("Cytotoxic", "NK 1"),
  "16" = c("naiveT", "naiveCD8 1"),
  "17" = c("naiveT", "naiveCD8 2"),
  "18" = c("naiveT", "naiveCD8 3"),
  "19" = c("naiveT", "naiveCD8 4"),
  "20" = c("Bcells", "Bcells 1"),
  "21" = c("Bcells", "Bcells 2"),
  "22" = c("Bcells", "Bcells 3"),
  "23" = c("Bcells", "plasmablasts"),
  "24" = c("Cytotoxic", "NK 2"),
  "25" = c("Cytotoxic", "NK 3"),
  "26" = c("other", "platelets 1"),
  "27" = c("other", "platelets 2"),
  "28" = c("other", "stressed"),
  "29" = c("other", "platelets 3"),
  "30" = c("Bcells", "Bcells 4"),
  "31" = c("myeloid", "classicalMono 1"),
  "32" = c("myeloid", "intermediateMono"),
  "33" = c("myeloid", "CD16p_Mono"),
  "34" = c("myeloid", "DCs"),
  "35" = c("myeloid", "classicalMono 2"),
  "36" = c("myeloid", "pDCs"),
  "37" = c("Bcells", "Bcells 5"),
  "38" = c("other", "progenitor cells")
) %>%
  pivot_longer(-granularity, names_to = "curated_clusters", values_to = "cluster_name") %>%
  pivot_wider(names_from = granularity, values_from = cluster_name) %>%
  mutate(curated_clusters = as.numeric(curated_clusters)) %>%
  mutate(cell_type = gsub(" [0-9]$", "", fine_clusters))


# build named palettes to ensure stability of cluster colors across all plots


#celltype
ct_names <- cluster_labels_table %>% pull(cell_type) %>% unique() %>% sort()

palette <- c(pal_rickandmorty()(12),pal_jco()(10), pal_futurama()(11), pal_tron()(7))

palette_ct <- palette[1:length(ct_names)]

names(palette_ct) <- ct_names

# coarse
palette_coarse <- pal_tron()(7) %>% darken(0.3)

names(palette_coarse) <- cluster_labels_table %>% pull(coarse_clusters) %>% unique() %>% sort()


# communities
set.seed(123)
palette_fine <- sample(rainbow(n = 38)) %>% desaturate(0.5)

names(palette_fine) <- cluster_labels_table %>% arrange(curated_clusters) %>% pull(fine_clusters)
palette_fine_n <- palette_fine
names(palette_fine_n) <- cluster_labels_table %>% arrange(curated_clusters) %>% pull(curated_clusters) %>% as.character()

palette_master <- c(palette_coarse, palette_ct, palette_fine, palette_fine_n)

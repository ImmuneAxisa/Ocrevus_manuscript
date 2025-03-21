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
library(ggridges)
library(here)

setwd(here())



#################################################
################ Import data ####################
#################################################

merged.seurat <- readRDS("results/PBMC/seurat_objects/master_seurat.rds")

merged.seurat <- AddMetaData(merged.seurat, 
                              PercentageFeatureSet(merged.seurat, pattern = "^MT-"), 
                              col.name = "percent.mt")


#################################################
######## filter out unwanted cells ##############
#################################################


# write a function to find the global maximum of the density
high_point <- function(x, high_only = F) {
  den <- density(x)
  high <- data.frame(x = den[["x"]],
                     y = den[["y"]]) %>%
    dplyr::filter(y == max(y)) %>%
    pull(x) %>%
    mean()
  
  if(!high_only) {
    if(mean(x) < high) {
      return(c(y = high,ymin = mean(x), ymax = high))
    } else {
      return(c(y = high,ymin = high, ymax = mean(x)))
    }
  } else {
    return(high)
  }
  
    
  
}


raw_meta <- merged.seurat@meta.data %>%
  tibble::rownames_to_column("cell.id") %>%
  separate(cell.id, c("patientID", "t" , "barcode"), sep = "_", remove = F) %>%
  mutate(t = ifelse(grepl("6m|6M", t), "t6M", "t0M")) %>%
  mutate(sample = paste(patientID, t, sep = "_"))
  
  
  
  
  
##---------------------------  
# plotting original data

p1 <- raw_meta %>%
  ggplot(aes(nCount_RNA, sample, fill = patientID)) +
  geom_density_ridges(quantile_lines = TRUE) +
  scale_x_log10(limits = c(NA,NA)) +
  stat_summary(fun.data = "high_point", geom = "pointrange", position = position_nudge(y = 0.8)) +
  scale_fill_cyclical(values = c("tan", "lightskyblue3")) +
  annotation_logticks(sides = "b") +
  geom_vline(xintercept = 1200) +
  ggtitle("All data, setting a lower bound \n to eliminate poor quality cells")

p2 <- raw_meta %>%
  dplyr::filter(nCount_RNA > 1200 & nFeature_RNA > 500 & nFeature_RNA < 4000) %>%
  ggplot(aes(nFeature_RNA, sample, fill = patientID)) +
  geom_density_ridges(quantile_lines = TRUE) +
  stat_summary(fun.data = "high_point", geom = "pointrange", position = position_nudge(y = 0.8)) +
  scale_x_log10() +
  scale_fill_cyclical(values = c("tan", "lightskyblue3")) +
  annotation_logticks(sides = "b") +
  ggtitle("After removal of bad quality cells and \n adding an upper bound to remove gross outliers")




##---------------------------  
# center with high point




test <- raw_meta %>%
  dplyr::filter(nCount_RNA > 1200 & nFeature_RNA > 500 & nFeature_RNA < 4000) %>%
  mutate(nCount_RNA = log10(nCount_RNA)) %>%
  mutate(nFeature_RNA = log10(nFeature_RNA)) %>%
  group_by(sample) %>%
  nest() %>%
  mutate(high_count = purrr::map(data, ~rep(high_point(.$nCount_RNA, high_only = T),
                                      nrow(.)))) %>%
  mutate(high_feature = purrr::map(data, ~rep(high_point(.$nFeature_RNA, high_only = T),
                                    nrow(.)))) %>%
  unnest() %>%
  ungroup() %>%
  mutate(nCount_RNA = 10^(nCount_RNA - high_count + mean(nCount_RNA))) %>%
  mutate(nFeature_RNA = 10^(nFeature_RNA - high_feature + mean(nFeature_RNA)))
  


# On counts

p3 <- test %>%
  ggplot(aes(nCount_RNA, sample, fill = patientID)) +
  geom_density_ridges(quantile_lines = TRUE) +
  stat_summary(fun.data = "high_point", geom = "pointrange", position = position_nudge(y = 0.8)) +
  scale_x_log10() +
  scale_fill_cyclical(values = c("tan", "lightskyblue3")) +
  annotation_logticks(sides = "b") +
  ggtitle("After removing outlier and \n centering on the global maxima")






# on features


p4 <- test %>%
  ggplot(aes(nFeature_RNA, sample, fill = patientID)) +
  geom_density_ridges(quantile_lines = TRUE) +
  stat_summary(fun.data = "high_point", geom = "pointrange", position = position_nudge(y = 0.8)) +
  scale_x_log10() +
  scale_fill_cyclical(values = c("tan", "lightskyblue3")) +
  annotation_logticks(sides = "b") +
  geom_vline(xintercept = c(900, 1800)) +
  ggtitle("After removing outlier and \n centering on the global maxima")
  
  
  
##---------------------------  
# apply filters and inspect output
  
# filter to only wanted cells
white_list <- test %>%
    dplyr::filter(nFeature_RNA > 900 & nFeature_RNA < 1800 & nCount_RNA > 2000 & nCount_RNA < 10000)

  
p5 <- white_list %>%
  ggplot(aes(nCount_RNA, sample, fill = patientID)) +
  geom_density_ridges(quantile_lines = TRUE) +
  stat_summary(fun.data = "high_point", geom = "pointrange", position = position_nudge(y = 0.8)) +
  scale_x_log10() +
  scale_fill_cyclical(values = c("tan", "lightskyblue3")) +
  annotation_logticks(sides = "b") +
  #geom_vline(xintercept = c(900, 1800)) +
  ggtitle("After removing tails with features cut offs")
  
  
p6 <- white_list %>%
  ggplot(aes(nFeature_RNA, sample, fill = patientID)) +
  geom_density_ridges(quantile_lines = TRUE) +
  stat_summary(fun.data = "high_point", geom = "pointrange", position = position_nudge(y = 0.8)) +
  scale_x_log10() +
  scale_fill_cyclical(values = c("tan", "lightskyblue3")) +
  annotation_logticks(sides = "b") +
  #geom_vline(xintercept = c(900, 1800)) +
  ggtitle("After removing tails with features cut offs")
  
  
##---------------------------   
# final plots: scatters

p7 <- white_list %>%
  ggplot(aes(nCount_RNA, nFeature_RNA)) +
  geom_hex(bins = 80, aes(color = ..count..)) +
  scale_x_log10() + 
  scale_y_log10() +
  annotation_logticks() +
  scale_fill_viridis() +
  scale_color_viridis()
 
  
p8 <- white_list %>%
  ggplot(aes(nCount_RNA, percent.mt)) +
  geom_hex(bins = 80, aes(color = ..count..)) +
  scale_x_log10() + 
  annotation_logticks(sides = "b") +
  scale_fill_viridis() +
  scale_color_viridis() +
  geom_hline(yintercept = 15, linetype = "dashed", color = "red")
  
  

##---------------------------  
# save plots

plot_list <- list(p1,p2,
                    p3,p4,
                    p5,p6,
                    p7,p8)

ggsave("results/PBMC/seurat_QC/master_seurat_QC_plots.pdf", 
        marrangeGrob(plot_list, nrow=1, ncol=2), 
        width = 30, height = 20, units = "cm")

#################################################
####### filter out the cells from seurat ########
#################################################

# lastly add a cut off for percentage of mitochondrial reads
white_list <- test %>%
    dplyr::filter(percent.mt < 15)


white_list <- pull(white_list, cell.id)


boolean <- row.names(merged.seurat@meta.data) %in% white_list


print("number of cells to keep / discard")
summary(boolean)

filtered.seurat <- merged.seurat[,boolean]

saveRDS(filtered.seurat, "results/PBMC/seurat_objects/master_filtered_seurat.rds")





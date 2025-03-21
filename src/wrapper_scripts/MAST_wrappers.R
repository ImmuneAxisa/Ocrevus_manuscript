############################################
############## MAST wrappers ###############
############################################




#################################################
################### Set up ######################
#################################################

## Take the input argument
# args <- commandArgs(TRUE)
# sample_name <- args[1]
# exp_matrix <- args[2]
# demux_best <- args[3]
# out_path <- args[4]

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
library(broom)
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
library(harmony)
library(future.apply)
library(speedglm)
library(rslurm)
library(here)
library(MAST)
library(qs)
library(ComplexHeatmap)


## ggplot theme

theme_mrl <- function(x = 1) {
  theme_minimal() +
    theme(
      axis.line = element_line(),
      axis.ticks.x = element_line(),
      axis.ticks.y = element_line(),
      axis.text.x = element_text(
        size = 12 * x,
        face = "bold",
        angle = 45,
        vjust = 0.6
      ),
      axis.text.y = element_text(size = 12 * x, face = "bold"),
      axis.title.x = element_text(size = 12 * x, face = "bold"),
      axis.title.y = element_text(size = 12 * x, face = "bold"),
      strip.background = element_rect(
        fill = "gray20",
        colour = "gray20",
        linetype = "solid"
      ),
      strip.text = element_text(
        size = 14 * x,
        colour = "white",
        face = "bold"
      ),
      legend.title = element_text(size = 14 * x, face = "bold"),
      legend.text = element_text(
        size = 12 * x,
        color = "gray20",
        face = "bold"
      ),
      legend.background = element_rect(fill = "transparent", colour = "transparent"),
      plot.title =  element_text(
        hjust = 0.5,
        vjust = 2,
        face = "bold"
      ),
      plot.subtitle = element_text(
        hjust = 0.5,
        vjust = 3,
        face = "italic"
      ),
      plot.caption = element_text(hjust = 0, face = "italic")
    )
}


catn <- function(...) {
  cat("\n", rep("#",50), "\n", sep = "")
  cat(...)
  cat("\n", rep("#",50), "\n", sep = "")
}


#################################################
################### functions ###################
#################################################

########################################################################
### cpc_function: return count per cell vector and print the distribution
### Return the cpc vector for filtering the GEX matrix
cpc_plot <- function(matrix, threshold = 0.01) {
  cpc <- Matrix::rowMeans(tmp_count)
  
  
  plot <- qplot(x = cpc + 10 ^ -4) +
    geom_vline(xintercept = threshold, color = "red") +
    labs(
      title = paste(
        "MAST on",
        clust,
        "with",
        ncol(tmp_count),
        "cells \n",
        paste(form, collapse = " ")
      ),
      caption = paste(sum(cpc > threshold), "genes pass of", length(cpc))
    ) +
    scale_x_log10()
  
  
  
  print(plot_grid(plot))
  
  
  return(cpc)
  
}

############################################################
### extract estimates and p values for the term of interest
get_MAST_results <- function(zlm_results, coeff) {
  # LRT test p value on reduced model
  pvalLRT <-
    lrTest(zlm_results, CoefficientHypothesis(name))[, , 3] %>% # 3rd dim contains pvalues
    as.data.frame() %>%
    tibble::rownames_to_column("gene_id") %>%
    rename_with( ~ paste0("pvalLRT_", .), .cols = -gene_id)
  
  # Wald test p value on coef
  pvalWald <-
    waldTest(zlm_results, CoefficientHypothesis(name))[, , 3] %>% # 3rd dim contains pvalues
    as.data.frame() %>%
    tibble::rownames_to_column("gene_id") %>%
    rename_with( ~ paste0("pvalWald_", .), .cols = -gene_id)
  
  # coeffs
  coeff_C <- coef(zlm_results, "C") %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene_id") %>%
    dplyr::select(gene_id, all_of(name)) %>%
    rename_with( ~ paste0("lfcC_", .), .cols = -gene_id)
  
  coeff_D <- coef(zlm_results, "D") %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene_id") %>%
    dplyr::select(gene_id, all_of(name)) %>%
    rename_with( ~ paste0("lfcD_", .), .cols = -gene_id)
  
  # SE of coeffs
  SEcoeff_C <- se.coef(zlm_results, "C") %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene_id") %>%
    dplyr::select(gene_id, all_of(name)) %>%
    rename_with( ~ paste0("SElfcC_", .), .cols = -gene_id)
  
  SEcoeff_D <- se.coef(zlm_results, "D") %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene_id") %>%
    dplyr::select(gene_id, all_of(name)) %>%
    rename_with( ~ paste0("SElfcD_", .), .cols = -gene_id)
  
  
  # assemble results into master table
  results <- pvalLRT %>%
    left_join(pvalWald) %>%
    left_join(coeff_C) %>%
    left_join(coeff_D) %>%
    left_join(SEcoeff_C) %>%
    left_join(SEcoeff_D) %>%
    left_join(data.frame("cpc" = cpc) %>%
                tibble::rownames_to_column("gene_id"))
  
  return(results)
  
}

#####################################################################
### creates some plots from the MAST results for the term of interest

plot_MAST_results <- function(results, name, partial_residuals, meta) {
  
  
  plot <- results %>%
    pivot_longer(-gene_id) %>%
    ggplot(aes(value)) +
    geom_histogram() +
    facet_wrap(. ~ name, scales = "free", ncol = 2) +
    theme_mrl()
  
  print(plot_grid(plot))
  
  catn("plot1")
  
  #continuous
  volcano_plot <- results %>%
    mutate(padj = p.adjust(pvalWald_hurdle)) %>%
    mutate(mlog10_p = -log10(pvalWald_hurdle + 10 ^ -250)) %>%
    ggplot(aes_string(paste0("lfcC_", name), "mlog10_p")) +
    geom_point(aes(color = padj < 0.1)) +
    geom_density2d() +
    theme_mrl()
  
  
  ma_plot <- results %>%
    mutate(padj = p.adjust(pvalWald_hurdle)) %>%
    mutate(mlog10_p = -log10(pvalWald_hurdle)) %>%
    ggplot(aes_string("cpc", paste0("lfcC_", name))) +
    geom_point(aes(color = padj < 0.1)) +
    geom_density2d() +
    scale_x_log10() +
    theme_mrl()
  
  print(plot_grid(volcano_plot, ma_plot, ncol = 1))
  
  catn("plot2")
  
  #discrete
  volcano_plot <- results %>%
    mutate(padj = p.adjust(pvalWald_hurdle)) %>%
    mutate(mlog10_p = -log10(pvalWald_hurdle + 10 ^ -250)) %>%
    ggplot(aes_string(paste0("lfcD_", name), "mlog10_p")) +
    geom_point(aes(color = padj < 0.1)) +
    geom_density2d() +
    theme_mrl()
  
  
  ma_plot <- results %>%
    mutate(padj = p.adjust(pvalWald_hurdle)) %>%
    mutate(mlog10_p = -log10(pvalWald_hurdle)) %>%
    ggplot(aes_string("cpc", paste0("lfcD_", name))) +
    geom_point(aes(color = padj < 0.1)) +
    geom_density2d() +
    scale_x_log10() +
    theme_mrl()
  
  print(plot_grid(volcano_plot, ma_plot, ncol = 1))
  
  catn("plot3")
  
  #heatmap
  degs <- results %>%
    mutate(padj = p.adjust(pvalWald_hurdle)) %>%
    dplyr::filter(
      padj < 0.1 &
        !grepl(
          "^RP[SL]|^MT-|^TRA[CVJ]|^TRB[CVDJ]|^TRD[CVJ]|^TRG[CVDJ]|^Gm[0-9]|Rik$",
          gene_id,
          ignore.case = T
        )
    )
  
  if(nrow(degs) > 0) {
    if (nrow(degs) > 1000) {
      degs <- degs %>%
        top_n(1000, -padj) %>%
        pull(gene_id)
      
    } else {
      degs <- degs %>%
        pull(gene_id)
    }
    
    meta_hm <- meta %>%
      dplyr::select(patient,
                    LP,
                    index) %>%
      tibble::column_to_rownames("index")
    
    km <- ceiling(length(degs)/300)
    
    catn("hm prep on", length(degs), "degs")
    
    
    
    hm <- partial_residuals[degs,] %>%
      t() %>%
      scale() %>%
      t() %>% {
        .[rowSums(is.na(.)) == 0,] # remove NA genes created by scaling low expression/no variance genes
      } %>% {
        ComplexHeatmap::Heatmap(
          .,
          #   right_annotation = rowAnnotation(link = anno_mark(
          #   at = grep(annotate, rownames(.)),
          #   labels = grep(annotate, rownames(.), value = T)
          # )),
          show_column_names = F,
          show_row_names = T,
          row_names_gp = gpar(fontsize = 6),
          top_annotation = HeatmapAnnotation(df = meta_hm[colnames(.),]),
          col = colorRamp2(quantile(., c(0.05, 0.5, 0.95), na.rm = T),
                           c("royalblue", "grey", "red4")),
          km = km
        )
      }
    
    print(hm)
    
    catn("hm")
    
  }
  
}

## ----library and source---------------------------------------------------------------------------
library(symphony)
library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(tidyverse)
library(Matrix)
library(sctransform)

setwd('/home/rstudio/autoimmune_10x')
source("./script/symphony/vignettes/utils_seurat.R")
source('./script/functions.R')

## ----parameter setting----------------------------------------------------------------------------
project.name <- "ocrevus_CSF_Yale"
prefix_data <- c("./data/")
prefix <- paste0("./output/", project.name, "/", project.name)

df.metadata <- read_csv(paste0(prefix_data,'/', project.name, '/', project.name, "_metadata.csv"))

## ----load h5seurat--------------------------------------------------------------------------------
# load(paste0(prefix_data, "ref_Reference_Mapping.RData"))
load(paste0("./data/ref_Reference_Mapping_20220525.RData"))

## ----process--------------------------------------------------------------------------------------
for (sample in df.metadata$sample) {
  print(sample)
  prefix = paste0("./output/", project.name, "/", sample)
  if (file.exists(paste0(prefix, "_CD4T_AssayData.rds"))) {
    query_obj <- readRDS(paste0(prefix, "_CD4T_AssayData.rds"))
    reference_mapping(ref, query_obj, prefix)
  } else {
    print("skipped...")
  }
}

df <- 0
for (sample in df.metadata$sample) {
  prefix = paste0("./output/", project.name, "/", sample)
  f <- paste0(prefix, "_queryL1_Reference_Mapping.csv")
  if (file.exists(f)){
    # print(f)
    d <- read.csv(f)
    if (df==0) {
      df <- d
    } else {
      df <- df %>% bind_rows(d)
    }
  }
}
write_csv(df, file = paste0("./output/", project.name, "/", project.name, "_queryL1_Reference_Mapping.csv"))

df <- 0
for (sample in df.metadata$sample) {
  prefix = paste0("./output/", project.name, "/", sample)
  f <- paste0(prefix, "_queryL2_Reference_Mapping.csv")
  if (file.exists(f)){
    # print(f)
    d <- read.csv(f)
    if (df==0) {
      df <- d
    } else {
      df <- df %>% bind_rows(d)
    }
  }
}
write_csv(df, file = paste0("./output/", project.name, "/", project.name, "_queryL2_Reference_Mapping.csv"))

df <- 0
for (sample in df.metadata$sample) {
  prefix = paste0("./output/", project.name, "/", sample)
  f <- paste0(prefix, "_Reference_Mapping.csv")
  if (file.exists(f)){
    # print(f)
    d <- read.csv(f)
    if (df==0) {
      df <- d
    } else {
      df <- df %>% bind_rows(d)
    }
  }
}
write_csv(df, file = paste0("./output/", project.name, "/", project.name, "_Reference_Mapping.csv"))

## ----save.image-----------------------------------------------------------------------------------
prefix = paste0("./output/", project.name, "/")
save.image(file = paste0(prefix, "220423_reference_mapping.RData"))

## ----make single h5ad-----------------------------------------------------------------------------------
vec.q <- c()
for (sample in df.metadata$sample) {
  print(sample)
  prefix_sample = paste0("./output/", project.name, "/", sample)
  if (file.exists(paste0(prefix_sample, "_CD4T_AssayData.rds"))) {
    count_obj <- readRDS(paste0(prefix_sample, "_CD4T_AssayData.rds"))
    meta_obj <- readRDS(paste0(prefix_sample, "_CD4T_MetaData.rds"))
    q <- CreateSeuratObject(counts=count_obj, assay="RNA", meta.data=meta_obj)
    vec.q <- append(vec.q, c(q))
    
  } else {
    print("skipped...")
  }
}
query_obj <- merge(vec.q[[1]], vec.q[2:length(vec.q)])
saveRDS(query_obj@meta.data, './output/ocrevus_CSF_Yale/ocrevus_CSF_Yale_CD4T_MetaData.rds')

SaveH5Seurat(query_obj, filename = './output/ocrevus_CSF_Yale/ocrevus_CSF_Yale_CD4T.h5Seurat')
Convert("./output/ocrevus_CSF_Yale/ocrevus_CSF_Yale_CD4T.h5Seurat", dest = "h5ad")


## ----sessionInfo----------------------------------------------------------------------------------
sessionInfo()

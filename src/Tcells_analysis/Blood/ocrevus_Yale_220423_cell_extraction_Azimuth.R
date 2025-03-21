## ----library and source---------------------------------------------------------------------------
library(Seurat)
library(SeuratData)
library(Azimuth)
library(patchwork)
library(tidyverse)
library(sctransform)
library(Matrix)
library(readxl)

setwd('/home/rstudio/autoimmune_10x')
source('./script/functions.R')


## ----parameter setting----------------------------------------------------------------------------
project.name <- c("ocrevus_Yale")
prefix_data <- paste0("./data/", project.name, "/")
prefix <- paste0("./output/", project.name, "/", project.name)

## ----Load_the_reference---------------------------------------------------------------------------
# reference <- LoadReference(path = "https://seurat.nygenome.org/azimuth/references/v1.0.0/human_pbmc")
reference <- LoadReference(path = "./data/azimuth/reference/v1.0.0/human_pbmc")

## ----reading data---------------------------------------------------------------------------------
query <- readRDS('./data/ocrevus_Yale/master_seurat_post_harmony_iterative_batch_emulsion.rds')
xlsx.meta <- read_xlsx('./data/ocrevus_Yale/Ocrevus 10x patient list_demographics.xlsx')
xlsx.meta <- xlsx.meta %>% rename("donor" = "MS ID", 
                                  "age" = "Age",
                                  "sex" = "Sex")
# merge(query@meta.data, xlsx.meta, by = "donor", all.x = TRUE)
query@meta.data <- query@meta.data %>% mutate(rownames = row.names(query@meta.data))
query@meta.data <- left_join(query@meta.data, xlsx.meta,  by="patient")
rownames(query@meta.data) <- query@meta.data$rownames

query@meta.data['disease'] = 'MS'
query@meta.data['disease_duration'] = NA
query@meta.data['condition'] <- recode(query@meta.data[['treatment']], 
                                       `tbase` = "baseline", `t0M` = "baseline", 
                                       `t06mo` = "post", `t6M` = "post")

query@meta.data['race'] = NA
query@meta.data['sex'] <- recode(query@meta.data[['sex']], 
                                       `F` = "female", `M` = "male")

query@meta.data['immunosuppressant'] <- recode(query@meta.data[['treatment']], 
                                               `tbase` = NA_character_, `t0M` = NA_character_, 
                                               `t06mo` = "ocrevus", `t6M` = "ocrevus")

query@meta.data['immunosuppressant_duration'] <- recode(query@meta.data[['treatment']], 
                                                        `tbase` = 0, `t0M` = 0, 
                                                        `t06mo` = 0.5, `t6M` = 0.5)

df.metadata <- query@meta.data
rownames(df.metadata) <- NULL
df.metadata <- df.metadata %>% select(c("donor", "treatment", "orig.ident", "sample", "batch", "disease", "disease_duration",
                                        "condition", "age", "race", "sex", "immunosuppressant", "immunosuppressant_duration")) %>%
  unique()
write_csv(df.metadata, paste0(prefix_data,'/', project.name, "_metadata.csv"))

## ----process--------------------------------------------------------------------------------------------------
for (sample in unique(query@meta.data$sample)) {
  print(sample)
  prefix = paste0("./output/", project.name, "/", sample)
  
  s = sample
  if (sum(query@meta.data['sample'] == s) > 0) {
    q <- subset(query, subset = sample == as.character(s))
    q <- CreateSeuratObject(counts = q@assays[["RNA"]]@counts,
                            project = project.name,
                            assay = "RNA",
                            meta.data = q@meta.data,
                            min.cells = 3,
                            min.features = 200) 
    extract_cells(q, reference, prefix)
  } else {
    print('No cell included in the sample. skipped...')
  }
}



## ----save.image-----------------------------------------------------------------------------------
prefix = paste0("./output/", project.name, "/")
save.image(file = paste0(prefix, "220423_cell_extraction_Azimuth.RData"))


## ----sessionInfo----------------------------------------------------------------------------------
sessionInfo()


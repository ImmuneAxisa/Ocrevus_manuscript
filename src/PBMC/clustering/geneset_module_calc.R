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
library(harmony)
library(presto)
library(here)



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
setwd(here())



objects_dir <- "results/PBMC/seurat_objects/"
dir.create(objects_dir, showWarnings = FALSE)


results_dir <- "results/PBMC/clustering"
dir.create(results_dir, showWarnings = FALSE)

#########################################################
###################### import data ######################
#########################################################


# get seurat object
merged.seurat <- readRDS(file.path(objects_dir,"curated_types_master_seurat_harmony_iterative_batch_emulsion.rds"))



# prepare genesets to test:
# TCR and IG gene segments
# cell types signatures derived from DICE online DGE portal

genesets <- list(

    TRAB = c("TRBV1","TRBV2","TRBV3-1","TRBV4-1","TRBV5-1","TRBV6-1","TRBV7-1","TRBV4-2","TRBV6-2","TRBV7-2","TRBV8-1","TRBV5-2","TRBV6-4","TRBV7-3","TRBV8-2","TRBV5-3","TRBV9","TRBV10-1","TRBV11-1","TRBV12-1","TRBV10-2","TRBV11-2","TRBV12-2","TRBV6-5","TRBV7-4","TRBV5-4","TRBV6-6","TRBV7-5","TRBV5-5","TRBV6-7","TRBV7-6","TRBV5-6","TRBV6-8","TRBV7-7","TRBV5-7","TRBV7-9","TRBV13","TRBV10-3","TRBV11-3","TRBV12-3","TRBV12-4","TRBV12-5","TRBV14","TRBV15","TRBV16","TRBV17","TRBV18","TRBV19","TRBV20-1","TRBV21-1","TRBV22-1","TRBV23-1","TRBV24-1","TRBV25-1","TRBVA","TRBV26","TRBVB","TRBV27","TRBV28","TRBV29-1","TRBJ1-1","TRBJ1-2","TRBJ1-3","TRBJ1-4","TRBJ1-5","TRBJ1-6","TRBC1","TRBJ2-1","TRBJ2-2","TRBJ2-2P","TRBJ2-3","TRBJ2-4","TRBJ2-5","TRBJ2-6","TRBJ2-7","TRBC2","TRBV30","TRBV20OR9-2","TRBV21OR9-2","TRBV22OR9-2","TRBV23OR9-2","TRBV24OR9-2","TRBV25OR9-2","TRBV26OR9-2","TRBV29OR9-2","TRAV1-1","TRAV1-2","TRAV2","TRAV3","TRAV4","TRAV5","TRAV6","TRAV7","TRAV8-1","TRAV9-1","TRAV10","TRAV11","TRAV12-1","TRAV8-2","TRAV8-3","TRAV13-1","TRAV12-2","TRAV8-4","TRAV8-5","TRAV13-2","TRAV14DV4","TRAV9-2","TRAV15","TRAV12-3","TRAV8-6","TRAV16","TRAV17","TRAV18","TRAV19","TRAV20","TRAV21","TRAV22","TRAV23DV6","TRAV24","TRAV25","TRAV26-1","TRAV8-7","TRAV27","TRAV28","TRAV29DV5","TRAV30","TRAV31","TRAV32","TRAV33","TRAV26-2","TRAV34","TRAV35","TRAV36DV7","TRAV37","TRAV38-1","TRAV38-2DV8","TRAV39","TRAV40","TRAV41","TRAJ61","TRAJ60","TRAJ59","TRAJ58","TRAJ57","TRAJ56","TRAJ55","TRAJ54","TRAJ53","TRAJ52","TRAJ51","TRAJ50","TRAJ49","TRAJ48","TRAJ47","TRAJ46","TRAJ45","TRAJ44","TRAJ43","TRAJ42","TRAJ41","TRAJ40","TRAJ39","TRAJ38","TRAJ37","TRAJ36","TRAJ35","TRAJ34","TRAJ33","TRAJ32","TRAJ31","TRAJ30","TRAJ29","TRAJ28","TRAJ27","TRAJ26","TRAJ25","TRAJ24","TRAJ23","TRAJ22","TRAJ21","TRAJ20","TRAJ19","TRAJ18","TRAJ17","TRAJ16","TRAJ14","TRAJ13","TRAJ12","TRAJ11","TRAJ10","TRAJ9","TRAJ8","TRAJ7","TRAJ6","TRAJ5","TRAJ4","TRAJ3","TRAJ2","TRAJ1","TRAC"),
    
    TRGD = c("TRGC2","TRGJ2","TRGJP2","TRGC1","TRGJ1","TRGJP","TRGJP1","TRGV11","TRGVB","TRGV10","TRGV9","TRGVA","TRGV8","TRGV7","TRGV6","TRG-AS1","TRGV5P","TRGV5","TRGV4","TRGV3","TRGV2","TRGV1","TRDMT1","TRDV1","TRDV2","TRDD1","TRDD2","TRDD3","TRDJ1","TRDJ4","TRDJ2","TRDJ3","TRDC","TRDV3"),
    
    IG = c("IGLV8OR8-1","IGHEP2","IGHMBP2","IGHA2","IGHE","IGHG4","IGHG2","IGHGP","IGHA1","IGHEP1","IGHG1","IGHG3","IGHD","IGHM","IGHJ6","IGHJ3P","IGHJ5","IGHJ4","IGHJ3","IGHJ2P","IGHJ2","IGHJ1","IGHD7-27","IGHJ1P","IGHD1-26","IGHD6-25","IGHD5-24","IGHD4-23","IGHD3-22","IGHD2-21","IGHD1-20","IGHD6-19","IGHD5-18","IGHD4-17","IGHD3-16","IGHD2-15","IGHD1-14","IGHD6-13","IGHD5-12","IGHD4-11","IGHD3-10","IGHD3-9","IGHD2-8","IGHD1-7","IGHD6-6","IGHD5-5","IGHD4-4","IGHD3-3","IGHD2-2","IGHD1-1","IGHV6-1","IGHVII-1-1","IGHV1-2","IGHVIII-2-1","IGHV1-3","IGHV4-4","IGHV7-4-1","IGHV2-5","IGHVIII-5-1","IGHVIII-5-2","IGHV3-6","IGHV3-7","IGHV3-64D","IGHV5-10-1","IGHV3-11","IGHVIII-11-1","IGHV1-12","IGHV3-13","IGHVIII-13-1","IGHV1-14","IGHV3-15","IGHVII-15-1","IGHV3-16","IGHVIII-16-1","IGHV1-17","IGHV1-18","IGHV3-19","IGHV3-20","IGHV3-21","IGHV3-22","IGHVII-22-1","IGHVIII-22-2","IGHV3-23","IGHV1-24","IGHV3-25","IGHVIII-25-1","IGHV2-26","IGHVIII-26-1","IGHVII-26-2","IGHV7-27","IGHV4-28","IGHVII-28-1","IGHV3-32","IGHV3-30","IGHVII-30-1","IGHV3-30-2","IGHV4-31","IGHVII-30-21","IGHV3-29","IGHV3-33","IGHVII-33-1","IGHV3-33-2","IGHV4-34","IGHV7-34-1","IGHV3-35","IGHV3-36","IGHV3-37","IGHV3-38","IGHVIII-38-1","IGHV4-39","IGHV7-40","IGHVII-40-1","IGHV3-41","IGHV3-42","IGHV3-43","IGHVII-43-1","IGHVIII-44","IGHVIV-44-1","IGHVII-44-2","IGHV1-45","IGHV1-46","IGHVII-46-1","IGHV3-47","IGHVIII-47-1","IGHV3-48","IGHV3-49","IGHVII-49-1","IGHV3-50","IGHV5-51","IGHVIII-51-1","IGHVII-51-2","IGHV3-52","IGHV3-53","IGHVII-53-1","IGHV3-54","IGHV4-55","IGHV7-56","IGHV3-57","IGHV1-58","IGHV4-59","IGHV3-60","IGHVII-60-1","IGHV4-61","IGHV3-62","IGHVII-62-1","IGHV3-63","IGHV3-64","IGHV3-65","IGHVII-65-1","IGHV3-66","IGHV1-67","IGHVII-67-1","IGHVIII-67-2","IGHVIII-67-3","IGHVIII-67-4","IGHV1-68","IGHV1-69","IGHV2-70D","IGHV3-69-1","IGHV1-69-2","IGHV1-69D","IGHV2-70","IGHV3-71","IGHV3-72","IGHV3-73","IGHV3-74","IGHVII-74-1","IGHV3-75","IGHV3-76","IGHVIII-76-1","IGHV5-78","IGHVII-78-1","IGHV3-79","IGHV4-80","IGHV7-81","IGHVIII-82","IGHV1OR15-9","IGHV1OR15-2","IGHV3OR15-7","IGHD5OR15-5A","IGHD4OR15-4A","IGHD3OR15-3A","IGHD2OR15-2A","IGHD1OR15-1A","IGHV1OR15-6","IGHD5OR15-5B","IGHD4OR15-4B","IGHD3OR15-3B","IGHD2OR15-2B","IGHD1OR15-1B","IGHV1OR15-1","IGHV1OR15-3","IGHV4OR15-8","IGHV1OR15-4","IGHV1OR16-1","IGHV1OR16-3","IGHV3OR16-9","IGHV2OR16-5","IGHV3OR16-15","IGHV3OR16-6","IGHV1OR16-2","IGHV3OR16-10","IGHV1OR16-4","IGHV3OR16-8","IGHV3OR16-12","IGHV3OR16-13","IGHV3OR16-11","IGHV3OR16-7","IGLJCOR18","IGLON5","IGHV1OR21-1","IGLVI-70","IGLV4-69","IGLVI-68","IGLV10-54","IGLV10-67","IGLVIV-66-1","IGLVV-66","IGLVIV-65","IGLVIV-64","IGLVI-63","IGLV1-62","IGLV8-61","IGLV4-60","IGLVIV-59","IGLVV-58","IGLV6-57","IGLVI-56","IGLV11-55","IGLVIV-53","IGLV5-52","IGLV1-51","IGLV1-50","IGLV9-49","IGLV5-48","IGLV1-47","IGLV7-46","IGLV5-45","IGLV1-44","IGLV7-43","IGLVI-42","IGLVVII-41-1","IGLV1-41","IGLV1-40","IGLVI-38","IGLV5-37","IGLV1-36","IGLV7-35","IGLV2-34","IGLV2-33","IGLV3-32","IGLV3-31","IGLV3-30","IGLV3-29","IGLV2-28","IGLV3-27","IGLV3-26","IGLVVI-25-1","IGLV3-25","IGLV3-24","IGLV2-23","IGLVVI-22-1","IGLV3-22","IGLV3-21","IGLVI-20","IGLV3-19","IGLV2-18","IGLV3-17","IGLV3-16","IGLV3-15","IGLV2-14","IGLV3-13","IGLV3-12","IGLV2-11","IGLV3-10","IGLV3-9","IGLV2-8","IGLV3-7","IGLV3-6","IGLV2-5","IGLV3-4","IGLV4-3","IGLV3-2","IGLV3-1","IGLJ1","IGLC1","IGLJ2","IGLC2","IGLJ3","IGLC3","IGLJ4","IGLC4","IGLJ5","IGLC5","IGLJ6","IGLC6","IGLJ7","IGLC7","IGLL1","IGLVIVOR22-1","IGLCOR22-1","IGLCOR22-2","IGLVIVOR22-2","IGLL5"),     
    
    Treg = c("LAYN","CCR3","HPGD","AC002331.1","FANK1","CCR8","LRRC32","TNFRSF9","AP001189.4","CEACAM4","FOXP3","C15orf53","RP11-62C7.2","CSF2RB","METTL7A","IKZF2","PTGIR","HNRNPA1P21","ST8SIA6","SEMA3G","AC017002.2","MKI67","ZNF365","ARAP3","ATP1B1","RASGRP4","TRIB1","AC017002.1","LGALS3","CNTNAP1","CD70","RTKN2","PI16","PTPLA","RP11-338N10.3","ADPRH","CCR5","CCR2","HLA-DQA2","HLA-DOA","DUSP4","MYBL2","MIR146A","GNG8","HLA-DQA1","IL12RB2","HLA-DQB1-AS1","C10orf128","WBSCR27","HLA-DRB5","MIR4435-1HG","HLA-DRB6","NCAPH","CCDC141","FUT7","JAKMIP1","HLA-DRB1","HLA-DQB1","AC021188.4","FCRL3","DFNB31","HLA-DQB2","SELP","RP11-603B24.1","IL17RB","VAV3","FCRL1","AC108463.2","E2F2","ITGAM","RP11-1399P15.1","IRF4","GPR55","PDE4A","LMNA","VWA7","FAM129A","CXCR6","XXYLT1-AS2","LINC00152","TPX2","PMAIP1","IL2RA","PRDM1","ZBTB32","SLC1A4","CIITA","SNORA14","RP11-390B4.3","HLA-DRA","PTMS","IL12A","RP11-536K7.5","NTNG2","HIP1","CD101","ARHGEF12","LGALS1","CTLA4","SUOX","ZNF532","C20orf197","IKZF4","STAM","IL1R1","DUSP10","TNFRSF18","HLA-DPA1","RP1-244F24.1","RP11-214O1.2","F5","NUSAP1","APOBEC3H","LAIR2","TIMD4","CCNG2","MGST2","HLA-DPB1","LAG3","GTSCR1","UTS2","TTN","C12orf75","NCR3","HLA-DMA","ITPRIPL1","SWAP70","NCF4","LTK","LPAL2","HLA-DMB","TK1","SPSB1","SHMT2","RP11-493L12.6","RP4-742C19.12","PHLDA1"),
    
    memCD4 = c("PTPN13","PDCD1","BAIAP3","CXCR3","NPDC1","TNFRSF4","RP11-403A21.2","PLEKHG3","MIAT","DUSP4","MYO1F","PLXND1","PRR5L","ADAM19","CXCR5","LINC00892","WDR86-AS1","KLRB1","CST7","ST8SIA1","CCR4","NCF4","SMCO4","ANK1","MAF","TIGIT","RP11-342D11.3","RP11-330A16.1","SYT11","GZMK","AC006369.2","NIPAL2","CYB561","HSF4","TBX21","FAM129A","CLCF1","UST","GPR15","RP5-1028K7.2","FBXL8","NOD2","IGFBP4","RP11-479G22.8","PREX1","RP11-367G6.3","SSBP4","CSF1","MCOLN2","DUSP5","ANXA2","LINC00963","THEM5","IGFBP3","PIEZO1","IFNG","STOM","MIR5587"),

    naiveCD4 = c("CD4","CD40LG","ANKRD55","EMR4P","CYSLTR1","IL2RA","DACT1","TSHZ2","RP11-223C24.1","CTSL","F5","GPR183","ADTRP","FRY","RP11-1399P15.1","RNF175","RTKN2","SNORA67_ENSG00000207419","AC020571.3","TIAM1","FHIT","STAP1","AC074289.1","ALS2CL","AP3M2","LPAR6","SORCS3","GLUL","EPHA4","GPA33","EMR1","ANK3","U8_ENSG00000238840","TSPAN18","CTLA4","TP53INP1","TNFSF10","PLCL1","SLAMF1","ZBTB7B","FAM13A","DTX1","TMEM45B"),
    
    naiveCD8 = c("RP11-291B21.2","CD8A","FBLN2","RP11-103C16.2","KLRK1","KLRC4-KLRK1","CD8B","RP11-277P12.20","RP11-277P12.9","SFRP5","CXCR3","S100B","DBN1","RP3-477M7.5","ANAPC1P1","CD248","PTK7","NT5E","SH3BGRL2","HLA-DOA","N4BP3","MIAT","YBX3","NPDC1","CRTAM","GZMK","SLC17A9","HMGN2P17","CTSW","REG4","PTPRK","NKG7","TTC24","RP11-535A5.1","SOAT2","CTC-501O10.1","CST7","NSG1","FGFBP2","CD79A","AC092580.4","AC093609.1","NRCAM","TSPAN32","MVB12B","MXRA8","MATK","NCR3","PRR5","ABCA2","AC116614.1","GSE1","TPST2","ZNF496","KCNN4","ASIC1"),
    
    mono1 = c("G0S2","IL1B","S100A9","S100A8","FCN1","VCAN","TNFAIP2","CD14","CD93","APOBEC3A","STAB1","PID1","MARCO","S100A12","AQP9","FAM198B","FOLR3","CCDC149","CSF3R","FPR1","CCR1","ASGR2","C15orf38","LYZ","CSTA","CD163","VSTM1","SIGLEC9","NME8","SERPINB2","FBN2","IL8","PLAUR","ANPEP","LILRA5","CD36","CACNA2D3","LIN7A","ZNF697","TYROBP","LGALS2","CTB-61M7.2","FCAR","NLRP12","GPR162","SERPING1","CPVL","CLEC12A","CDA","PYGL","THBD","NRG1","SERPINA1","HBEGF","CD300E","CLEC7A","HK3","RAB32","CES1","PROK2","SIGLEC22P","RHOU","TMEM150B","SIGLEC1","SLCO3A1","CUEDC1","RPH3A","CST3"),
    
    mono2= c("CLEC4F","FCGR3A","RP11-290F20.3","AC104809.3","MEG3","AC104809.4","RP11-362F19.1","SH2D1B","SERPINA1","CYP4F22","LILRA5","C1QA","LYPD2","TMTC1","C1QB","CEACAM3","FCER1G","CLEC12A","HK3","EPB41L3","CUEDC1","RP11-1008C21.1","FCGR3B","ABCC3","CTSL","MAFB","SLC8A1"),
    
    NK = c("COL13A1","CCNJL","SH2D1B","LINGO2","GRIK4","ZMAT4","FAT4","BNC2","PRICKLE2","PCDH1","RNF165","LIM2","LINC00298","KRT86","KIR2DL1","NCR1","GZMB","CTD-3187F8.14","NCAM1","SIGLEC7","KIR2DP1","KIR2DS4","C19orf84","TRDJ1","RP11-1055B8.4","LINC00299","KIR3DL1","KIR2DL4","CMKLR1","CTD-2616J11.14","KRT81","MGAM","ADAMTS1","RP11-438B23.2","NMUR1","CXCR2","PRSS23","CXCR1","SIGLEC19P","KIFC3","CHSY3","FCGR3A","CTD-2616J11.3","KIR2DL3","ST8SIA6","RP11-845M18.6","CTD-2616J11.2","NIFKP6","TRDD3","AC007278.3","ATP8B4","DAB2","PRSS57","MIR4772","B3GNT7","RP11-705C15.5","IL18RAP","S1PR5","SLCO4C1","CX3CR1","GNLY","KIR3DX1","GPR56","SLC1A7","FASLG","FABP6","FHAD1","LINC00504","PODN","C2orf48","AC007278.2","CTD-2154I11.2","KLRF1","FCGR3B","NME8","SMPD4P1","CATSPER1","VSTM4","RP11-627G23.1","NCAM1-AS1","B3GAT1","RAMP1","GLI3","MLC1","BFSP1","PTPRB","IGFBP7","OLA1P2","U91324.1","GPR141","ATP9A","LGALS9B","RP13-39P12.3","EFNA5","RP11-1220K2.2")

)

#####################################
####### run module enrichments ######
#####################################

merged.seurat <- AddModuleScore(merged.seurat, 
                                genesets,
                                name = paste0(names(genesets),"xxx"))
                                
                                
new_cols <- sub("xxx.*","",colnames(merged.seurat@meta.data))

colnames(merged.seurat@meta.data) <- new_cols


meta <- merged.seurat@meta.data %>%
            tibble::rownames_to_column("cell.id")

saveRDS(meta, file.path(results_dir,"master_seurat_curated_types_table_harmony_iterative_batch_emulsion.rds"))

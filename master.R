##############################################################################
##############################################################################
############################## project master ################################
##############################################################################
##############################################################################

# Warning:
# This scripts runs through all project analyses, but is not executed as is.
# It is not intended as a stand-alone, unsupervised piece of code, but a documentation of the various procedures

library(here)
library(callr)
library(rslurm)


# Using HPC resources to run the analysis:
# ----------------------------------------
# Parts of the analysis are preferably run on an HPC cluster using a slurm scheduler
# we use a convenience function here that lets the user decide whether to
# run the script in the current session (default) or submit a job using the {rslurm} package. 
# If doing so you will need adapt slurm options to your HPC infrastructure.
# see rslurm::slurm_call function details on how to manage this.

rscript2 <- function(script_path, sbatch = FALSE, ...){
  
  if(sbatch) {
    slurm_call(callr::rscript, params = list(script = script_path), ...)
    message("job submitted to slurm scheduler, wait for its completion before running other analyses that depends on this script's results")
  } else {
    rscript(script_path)
  }
  return(invisible(NULL))
}

##############################################################################
##############################################################################
########################### Main PBMC analysis ###############################
##############################################################################
##############################################################################

####################################
########## preprocessing ###########
####################################

# create seurat object from cellranger outputs
rscript2("src/PBMC/GEX_QC/create_master_seurat.R")

# QC: remove poor quality and doublet cells
rscript2("src/PBMC/GEX_QC/QC_master_seurat.R")


####################################
##### Reductions and clustering ####
####################################

# PCA/harmony
# -----------
rscript2("src/PBMC/PCA/harmony_iterative_batch_emulsion.R", TRUE,
         slurm_options = list(nodes = 1, 
                              cpus_per_node = 1, 
                              mem = "119gb", 
                              partition = "general", 
                              error = here("slurm_%a.err"),
                              time = "12:00:00"),
         submit = T)

# Clustering and UMAPs/tSNE
# -------------------------

# First pass analysis at low resolution
rscript2("src/PBMC/clustering/clustering_iterations.R", TRUE,
         slurm_options = list(nodes = 1, 
                              cpus_per_node = 1, 
                              mem = "64gb", 
                              partition = "general", 
                              error = here("slurm_%a.err"),
                              time = "36:00:00"),
         submit = T)

# subclustering
rscript2("src/PBMC/clustering/subclustering_iterations.R", TRUE,
         slurm_options = list(nodes = 1, 
                              cpus_per_node = 1, 
                              mem = "64gb", 
                              partition = "general", 
                              error = here("slurm_%a.err"),
                              time = "36:00:00"),
         submit = T)



# MELD/PHATE
# ----------
# this is admittedly a bit convoluted, but hey it works!
# phate_meld.R prepares PC loadings and metadata tables to be given as input 
# for MELD and PHATE. In turn, it calls a wrapper shell script
# "src/wrapper_scripts/phate_meld_sbatcher_w_args.sh"
# which prepares slurm scripts to be submitted to the scheduler.

# Depending on your shell config, the shell called within R might not
# be able to call sbatch, in which case you will need 
# to submit the slurm script "manually" on your HPC cluster.
# If you do not wish to use slurm scheduler for those analyses, 
# the generated shell scripts can simply be executed.

# The submitted scripts call "src/wrapper_scripts/MELD_v2.1.py" and
# "src/wrapper_scripts/PHATE_v2.py" to run the analyses.

# to run those analyses you will need to install PHATE and MELD
# python packages. In our analyses this is done by setting up a
# virtual environment "MELD_env" through miniconda.
# the configuration of the conda environment can be found here:
# "src/wrapper_scripts/MELD_env_requirements.txt"

rscript2("src/PBMC/phate_meld/phate_meld.R")


####################################
###### biological annotations ######
####################################

# we assign cell types to clusters
# we computed some gene signatures and calls to reference
# to assign broad cell types

# clustering curation
# -------------------

# basic curation
rscript2("src/PBMC/clustering/clustering_manual_cleaning.R", TRUE,
         slurm_options = list(nodes = 1, 
                              cpus_per_node = 1, 
                              mem = "64gb", 
                              partition = "general", 
                              error = here("slurm_%a.err"),
                              time = "06:00:00"),
         submit = T)

# multi-level curation with coarse vs fine-grained types
# this is largely redundant with "clustering_manual_cleaning.R"
# but the 2 are kept because they generate slightly different outputs
# which are used subsequently in the analysis.
rscript2("src/PBMC/clustering/clustering_granularity_cleaning.R", TRUE,
         slurm_options = list(nodes = 1, 
                              cpus_per_node = 1, 
                              mem = "64gb", 
                              partition = "general", 
                              error = here("slurm_%a.err"),
                              time = "06:00:00"),
         submit = T)


# frequency variations with sccomp
# lotsa clusters, lotsa runtime
rscript2("src/PBMC/clustering/freq_variation_sccomp.R", TRUE,
         slurm_options = list(nodes = 1, 
                              cpus_per_node = 1, 
                              mem = "64gb", 
                              partition = "general", 
                              error = here("slurm_%a.err"),
                              time = "12:00:00"),
         submit = T)


# cluster markers 
# ---------------
rscript2("src/PBMC/clustering/cluster_markers.R", TRUE,
         slurm_options = list(nodes = 1, 
                              cpus_per_node = 1, 
                              mem = "64gb", 
                              partition = "general", 
                              error = here("slurm_%a.err"),
                              time = "06:00:00"),
         submit = T)

# geneset signatures and reference calls
# --------------------------------------
# custom geneset scores
rscript2("src/PBMC/clustering/geneset_module_calc.R")

# singleR on clusters
rscript2("src/PBMC/clustering/cell_typing_multi_res.R", TRUE,
         slurm_options = list(nodes = 1, 
                              cpus_per_node = 1, 
                              mem = "64gb", 
                              partition = "general", 
                              error = here("slurm_%a.err"),
                              time = "06:00:00"),
         submit = T)

# singleR at single cell level
rscript2("src/PBMC/clustering/cell_typing_sc.R", TRUE,
         slurm_options = list(nodes = 1, 
                              cpus_per_node = 1, 
                              mem = "64gb", 
                              partition = "general", 
                              error = here("slurm_%a.err"),
                              time = "24:00:00"),
         submit = T)


####################################
### Differential gene expression ###
####################################

# NEBULA run for all clusters
rscript2("src/PBMC/DGE/nebula_DGE_curated_clusters.R", TRUE,
         slurm_options = list(nodes = 1, 
                              cpus_per_node = 1, 
                              mem = "119gb", 
                              partition = "general", 
                              error = here("slurm_%a.err"),
                              time = "72:00:00"),
         submit = T)

# Second pass on CD16 monocytes (cluster33) to remove 
# contaminant cells polluting the DGE results
rscript2("src/PBMC/DGE/nebula_DGE_cluster33_noDoublets.R", TRUE,
         slurm_options = list(nodes = 1, 
                              cpus_per_node = 1, 
                              mem = "119gb", 
                              partition = "general", 
                              error = here("slurm_%a.err"),
                              time = "12:00:00"),
         submit = T)


####################################
# Results visualization & figures ##
####################################

# first pass analysis and QC
rmarkdown::render("doc/PBMC_initial_QC_and_exploration.Rmd")

# Advanced results and manuscript figure building
rmarkdown::render("doc/PBMC_analysis.Rmd")

# DGE and GSEA results and manuscript figure building
rmarkdown::render("doc/PBMC_DGE_GSEA.Rmd")

##############################################################################
##############################################################################
########################### Main CSF analysis ################################
##############################################################################
##############################################################################

# since the CSF dataset is relatively small, most analyses are done interactively in an Rmd file

####################################
####### Interactive analysis #######
####################################

# QC and basic analysis: dimred, clustering
rmarkdown::render("doc/CSF_MSandHA_Integrated_QC.Rmd")

# CSF signature in Macrophage and T cells
# refine to separate mixed lineage cluster
rmarkdown::render("doc/CSF_macrophage_cluster_inspection.Rmd")

# Macrophage exploration
rmarkdown::render("doc/CSF macrophage phenotype.Rmd")

# DC exploration
rmarkdown::render("doc/CSF DC pDC.Rmd")

####################################
##### Reductions and clustering ####
####################################

# cluster markers
rscript2("src/CSF_analysis/clustering_typing/cluster_markers.R")

# UMAPs
rscript2("src/CSF_analysis/clustering_typing/CSF_UMAPs.R")

# Phate/MELD: see similar PBMC analysis 
# above for explanation on how this is run
rscript2("src/CSF_analysis/phate_meld/phate_meld_CSF.R")

# freq variation with sccomp
rscript2("src/CSF_analysis/clustering_typing/freq_variation_sccomp_CSF.R")

####################################
### Differential gene expression ###
####################################

# DGE with DESeq2
rscript2("src/CSF_analysis/DGE/DESEq2_glmGamPois.R", TRUE,
         slurm_options = list(nodes = 1, 
                              cpus_per_node = 1, 
                              mem = "64gb", 
                              partition = "general", 
                              error = here("slurm_%a.err"),
                              time = "05:00:00"),
         submit = T)


# plot DGE results and GSEA against PBMC results 
rmarkdown::render("doc/CSF_DGE_GSEA.Rmd")

####################################
########### Paper figures ##########
####################################

rmarkdown::render("doc/figures_CSF.Rmd")


##############################################################################
##############################################################################
############################ T cells analysis ################################
##############################################################################
##############################################################################

# T cells analysis were run CSF and PBMC.
# releveant scripts are located here: src/Tcells_analysis
# may need some some tweaking of relative path to run.


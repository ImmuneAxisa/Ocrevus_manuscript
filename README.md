# Ocrevus_manuscript

This repository contains the code associated with the manuscript [Transcriptomic profiling after B-cell depletion reveals central and peripheral immune cell changes in multiple sclerosis](https://doi.org/10.1172/JCI182790). Analyzed Seurat objects are available on [figshare](https://doi.org/10.6084/m9.figshare.28204532

# Disclaimer

* Throughout this project, We use the `SLURM` job scheduling system to run some of the "heavy" analyses on the Yale University HPC clusters.
* The code and additional resources on figshare.com are provided as-is. We have validated that this runs without errors on the HPC cluster and on local MacOS systems, but do not guarantee that this will be the case in different computing environments.
* Main analyses portion of this code are written with relative paths that have as root project this directory, but some path might need to be adjusted for the T cell analysis.


# Content

* `master.R` details all procedures run throughout the project to generate the findings presented in the manuscript. It calls scripts stored in `src` and `doc` and is used to document the analysis steps.
* `src` contains all analyses steps
* `doc` contains `Rmd` documents used for data exploration and visualization.
* `prepare_final_objects.R` details the procedure to generate the objects uploaded on figshare.
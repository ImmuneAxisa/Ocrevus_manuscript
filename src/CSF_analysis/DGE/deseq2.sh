#!/bin/bash
#SBATCH --partition=day
#SBATCH --job-name=deseq
#SBATCH --time=5:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64gb
#SBATCH -o /gpfs/gibbs/project/hafler/pa326/Ocrevus/logs/%x.out
#SBATCH -e /gpfs/gibbs/project/hafler/pa326/Ocrevus/logs/%x.err


module load R/4.3.2-foss-2022b
module load R-bundle-Bioconductor/3.18-foss-2022b-R-4.3.2 
module load R-bundle-CRAN/2023.12-foss-2022b

Rscript /gpfs/gibbs/project/hafler/pa326/Ocrevus/src/CSF_analysis/DGE/DESEq2_glmGamPois.R
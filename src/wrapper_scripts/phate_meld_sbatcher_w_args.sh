#!/bin/bash

###############################################################
######################### set up ##############################
###############################################################


# command line arguments
# column names for the metadata table to be used in MELD
base_dir=$1 # project dir
results_dir=$2 # should contain the "phate_meld_input_csv" DIR to store csv inputs
input_csv_file=$3      # input files basenames (has to finish in _pca.csv and _meta.csv)
sample_labels=$4       # eg "treatment"
cluster=$5            # eg "cluster_label"
out_suffix=$6 # concat to sample name



# Raw data and working directories:
data_dir="${base_dir}/data"
input_dir="${results_dir}/phate_meld_input_csv"

# Output directories:
slurm_dir="${base_dir}/slurm"
logs_dir="${base_dir}/logs"
doc_dir="${base_dir}/doc"
results_dir="${base_dir}/phate_meld_output"


# Create any missing directories
directories=( $data_dir $logs_dir $doc_dir $slurm_dir \
  $input_dir $results_dir)
for directory in ${directories[@]}; do
  if [ ! -d $directory ]; then
    echo "Creating directory ${directory}"
    mkdir -p $directory
  fi
done

###############################################################
######################### sbatch ##############################
###############################################################

# Sample table prefixes: 
# structure must be <prefix>_pca.csv and <prefix>_meta.csv
# for both files the first column must be the cell IDs eg barcodes

# Array of samples with prefixes to be processed
#samples=(CD4 CD8 CD8NK B M)
samples=$input_csv_file
#samples=(B)




for sample in ${samples[@]} #array is deprecated but still functional with one argument
do
    # format inputs and check they exists
    pca=${input_dir}/${sample}_pca.csv
    meta=${input_dir}/${sample}_meta.csv
    
    if [ ! -f "$pca" ]; then
        echo "$pca not found"
        exit
    fi
    
    if [ ! -f "$meta" ]; then
        echo "$meta not found"
        exit
    fi
    
    sample=${sample}_${out_suffix}
    
    echo "processing $sample"
    
    #slurm scripts
    phate_script=${slurm_dir}/${sample}_phate.slurm
    meld_vfc_script=${slurm_dir}/${sample}_meld_vfc.slurm
  
  
############ PHATE ##############

    echo \
"#!/bin/bash
#SBATCH --partition=general
#SBATCH --job-name=${sample}.phate
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=119gb
#SBATCH --time=48:00:00
#SBATCH -o ${logs_dir}/${sample}.phate.out
#SBATCH -e ${logs_dir}/${sample}.phate.err

module load miniconda

source activate MELD_env

python -u ${base_dir}/src/wrapper_scripts/PHATE_v2.py \\
    -i ${pca} \\
    -o ${results_dir} \\
    -p ${sample}


source deactivate

sstat \${SLURM_JOBID}.batch


" > $phate_script

sbatch $phate_script


############ MELD ##############

    echo \
"#!/bin/bash
#SBATCH --partition=bigmem
#SBATCH --job-name=${sample}.meld_vfc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=500gb
#SBATCH --time=48:00:00
#SBATCH -o ${logs_dir}/${sample}.meld_vfc.out
#SBATCH -e ${logs_dir}/${sample}.meld_vfc.err


module load miniconda

conda init

source activate MELD_env

python -u ${base_dir}/src/wrapper_scripts/MELD_v2.1.py \\
    -i ${pca} \\
    -m ${meta} \\
    -l ${sample_labels} \\
    -c ${cluster} \\
    -o ${results_dir} \\
    -p ${sample}

source deactivate

sstat \${SLURM_JOBID}.batch

" > $meld_vfc_script

sbatch $meld_vfc_script


done
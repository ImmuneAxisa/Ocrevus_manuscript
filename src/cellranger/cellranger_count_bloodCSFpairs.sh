#!/bin/bash

###############################################################
######################### set up ##############################
###############################################################

# Raw data directories

raw_data="/gpfs/ycga/scratch60/hafler/pa326/Ocrevus/jenna/rawdata_HA"

#Reference genome
HHT="/gpfs/ycga/datasets/genomes/10xgenomics/refdata-cellranger-GRCh38-3.0.0"

#working directories:
base_dir="/home/pa326/project/Ocrevus"
data_dir="${base_dir}/data"
scratch_dir="/gpfs/ycga/scratch60/hafler/pa326/Ocrevus/cellranger"

# Output directories:
slurm_dir="${base_dir}/slurm"
logs_dir="${base_dir}/logs"
doc_dir="${base_dir}/doc"
results_dir="${base_dir}/results/"


# Create any missing directories
directories=( $data_dir $logs_dir $doc_dir $slurm_dir \
  $results_dir $scratch_dir)
for directory in ${directories[@]}; do
  if [ ! -d $directory ]; then
    echo "Creating directory ${directory}"
    mkdir -p $directory
  fi
done


###############################################################
#################### cellranger jobs ##########################
###############################################################

# Get list of unique stems to process

if ls -1 $raw_data | grep -q fastq; then
    
    echo "processing as all fastqs in the same directory"
    ls -1 $raw_data | grep HHT | grep fastq | sed 's/_S[0-9]*_L[0-9]*_[IR][1-2]_[0-9]*\..*//' | sort | uniq > ${scratch_dir}/sample.lst
    
else
    
    echo "processing as each sample fastq in separate directories"
    
    ls -1 $raw_data | sort | uniq > ${scratch_dir}/sample.lst
    
fi


# Loop through each sample to submit slurm jobs

while read sample

do 
    
    # if fastq in separate directories, append the directory name for each sample
    # else just use the path where all fastqs are stored
    if  ls -1 $raw_data | grep -q fastq; then
        raw_path=$raw_data
    else
        raw_path=${raw_data}/${sample}
        sample=$(ls -1 $raw_path | grep fastq | sed 's/_S[0-9]*_L[0-9]*_[IR][1-2]_[0-9]*\..*//' | sort | uniq)
        
        if [[ $sample == *" "* ]]; then
                echo "several samples in the same folder despite one sample per folder guessed for $sample"
                exit 1
        fi
    fi
    
    
    
    echo processing $sample
    
    # sample name for output directory
    sample_name=$(echo $sample | sed 's/_HHT|_MMT|_5P//')
    
    # slurm script file
    slurm_script=${slurm_dir}/${sample}.slurm
    
    
    echo \
"#!/bin/bash
#SBATCH --partition=general
#SBATCH --job-name=${sample}.cellranger
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=119gb
#SBATCH --time=24:00:00
#SBATCH -o ${logs_dir}/${sample}.cellranger.out
#SBATCH -e ${logs_dir}/${sample}.cellranger.err


# load cellranger
module load cellranger/3.1.0

cd ${scratch_dir}

cellranger count --id=${sample_name} \\
                        --transcriptome=$HHT \\
                        --fastqs=$raw_path \\
                        --sample=$sample

" > $slurm_script

    # submit slurm script directly if sure it'll work
    sbatch $slurm_script

done < ${scratch_dir}/sample.lst



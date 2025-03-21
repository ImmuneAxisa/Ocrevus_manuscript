#!/bin/bash

###############################################################
######################### set up ##############################
###############################################################

# Raw data directories
#raw_data="/ycga-gpfs/project/fas/lsprog/tools/external/data/PqR66O9Y5nmRYMtXzvruqvK2HZDoY/042221/"
#raw_data="/ycga-gpfs/project/fas/lsprog/tools/external/data/30cBhe3mdn06LcvORZcBr2Ys2fvxs/012621/"
raw_data="/ycga-gpfs/project/fas/lsprog/tools/external/data/IJKOfwvIeKDjr6xNSDzyjbJCVigh9/033021/"
#raw_data="/ycga-gpfs/project/fas/lsprog/tools/external/data/yWu101tGqtZnZPetMsAIxNkePvDsd/062121/"

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
    ls -1 $raw_data | grep HHT | grep fastq | sed 's/HHT.*/HHT/' | sort | uniq > ${scratch_dir}/sample.lst
    
else
    
    echo "processing as each sample fastq in separate directories"
    
    ls -1 $raw_data | grep HHT | grep fqs | sort | uniq > ${scratch_dir}/sample.lst
    
fi


# Loop through each sample to submit slurm jobs

while read sample

do 
    
    # if fastq in separate directories, append the directory name for each sample
    # else just use the path where all fastqs are stored
    if  [ -d "${sample}" ] ; then
        raw_path=${raw_data}/${sample}
        
    else
        raw_path=$raw_data
    fi
    
    sample=$(echo $sample | sed 's/_HHT_fqs//')
    
    echo processing $sample
    
    # sample name for output directory
    sample_name=$(echo $sample | sed 's/_HHT//')
    
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
    #sbatch $slurm_script

done < ${scratch_dir}/sample.lst



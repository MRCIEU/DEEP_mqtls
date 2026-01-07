#!/bin/bash -l
#SBATCH --job-name=deep_01c_ewas
#SBATCH --output=../job_reports/deep_01c_ewas_%j
#SBATCH --partition=compute 
#SBATCH --account=SSCM029144
#SBATCH --mem=128GB
#SBATCH --cpus-per-task=16
#SBATCH --time=36:00:00

# if you are in DEEP_mqtls/submission_scripts_template folder
cd ..
module load git

./01c-check_phenotypes_and_methylation.sh -c /user/work/er20212/alspac_QC_phase1/config_alspac

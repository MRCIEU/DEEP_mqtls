#!/bin/bash -l
#SBATCH --job-name=deep_01c_ewas
#SBATCH --output=../job_reports/deep_01c_aggregate_%j
#SBATCH --partition gpu,cpu
#SBATCH --mem=64GB
#SBATCH --ntasks=8
#SBATCH --time=6:0:0

# if you are in DEEP_mqtls/submission_scripts_template folder
cd ..
module load git

./01c-check_phenotypes_and_methylation.sh -c /path/to/config/file ewas

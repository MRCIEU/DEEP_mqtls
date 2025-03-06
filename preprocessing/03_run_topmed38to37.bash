#!/bin/bash
#SBATCH --job-name=38to37
#SBATCH --partition=mrcieu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=36:00:00
#SBATCH --mem=90G
#SBATCH --account=SSCM029144
#SBATCH --array=1-23%1
#SBATCH --output=job_reports/slurmTopMed-%A_%a.out

module load languages/R/4.4.1

INPUT_DIR="/user/work/er20212/data/"
INPUT_FILE="temp_${SLURM_ARRAY_TASK_ID}.snplist"

Rscript topmed38to37.R "${INPUT_DIR}${INPUT_FILE}"

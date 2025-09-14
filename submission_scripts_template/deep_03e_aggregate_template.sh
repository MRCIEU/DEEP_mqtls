#!/bin/bash -l
#SBATCH --job-name=deep_03e_aggregate
#SBATCH --output=../job_reports/deep_03e_aggregate_%j
#SBATCH --partition gpu,cpu
#SBATCH --mem=64GB
#SBATCH --ntasks=8
#SBATCH --time=6:0:0

cd ..
source resources/setup.sh "$@"
set -- $concatenated

bash ./resources/methylation/aggregate_adjustment2.sh

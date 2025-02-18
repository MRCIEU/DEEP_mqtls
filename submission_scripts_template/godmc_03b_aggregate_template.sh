#!/bin/bash -l
#SBATCH --job-name=GoDMC_03b_aggregate
#SBATCH --output=../job_reports/GoDMC_03b_aggregate_%j
#SBATCH --partition gpu,cpu
#SBATCH --mem=64GB
#SBATCH --ntasks=8
#SBATCH --time=6:0:0

cd ..
source resources/setup.sh "$@"
set -- $concatenated
bash ./resources/methylation/aggregate_adjustment1.sh

count0=`ls ${transformed_methylation_adjusted}* | grep ${transformed_methylation_adjusted}.Female.chrX | wc -l`
count1=`ls ${transformed_methylation_adjusted}* | grep ${transformed_methylation_adjusted}.Male.chrX | wc -l`
count2=`ls ${transformed_methylation_adjusted}* | grep ${transformed_methylation_adjusted}.Male.chrY |wc -l`

if [ ${count0} -gt 0 ]
then
    bash ./resources/methylation/aggregate_adjustment1.sh Female.chrX
else
    echo "No chrX CpGs from female samples"
fi

if [ ${count1} -gt 0 ]
then
    bash ./resources/methylation/aggregate_adjustment1.sh Male.chrX
else
    echo "No chrX CpGs from male samples"
fi

if [ ${count2} -gt 0 ]
then
    bash ./resources/methylation/aggregate_adjustment1.sh Male.chrY
else
    echo "No chrY CpGs from male samples"
fi

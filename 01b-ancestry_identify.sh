#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

exec &> >(tee ${section_01c_logfile})
print_version

# calculate maf from bfile
echo "Calculating MAF from bfile"
plink2 --bfile "${bfile_raw}" --freq --out "${section_01_dir}/"

#echo "Converting PLINK files to VCF format"
#plink2 --bfile "${bfile_raw}" \ 
#	   --export vcf \
#       --out "${vcf_file}"

echo "Running global PCA"
# global pca plot
${Python_directory}python "${scripts_directory}/resources/datacheck/ancestry_infer.py" \
    "${section_01_dir}/logs_b/hail.log" \
    "${bfile_raw}" \
	"${genome_build}" \
    "${study_name}" \
    "${scripts_directory}"

# local pca plot to be done
echo "Successfully completed script 1b"

#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

exec &> >(tee ${section_01c_logfile})
print_version

echo "Converting PLINK files to VCF format"
plink2 --bfile "${bfile_raw}" \
	   --export vcf \
       --out "${vcf_file}"

echo "Running global PCA"
# global pca plot
python "${scripts_directory}/resources/datacheck/ancestry_infer.py" \
    "${section_01_dir}/logs_b/hail.log" \
    "${vcf_file}.vcf" \
	"${genome_build}" # Use the genome build variable from the config file

# local pca plot
# to be done
echo "Successfully completed script 1b"
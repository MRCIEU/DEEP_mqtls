#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

exec &> >(tee ${section_01c_logfile})
print_version

echo "Converting PLINK files to VCF format"
plink2 --bfile "${bfile_raw}" \
	   --export vcf \
       --out "${vcf_file}"

# # Check genome build and liftover to 38 if build is 37
# ${R_directory}Rscript resources/datacheck/liftover.R \
# 	${bfile_raw} \
# 	${genome_build} \
# 	${miss_liftover}

# if [ -f ${miss_liftover} ]; then
#     plink2 --bfile "${bfile_raw}" \
#            --exclude-snps ${miss_liftover} \
#            --make-bed \
#            --out {liftover_data}

# 	plink2 --bfile "${bfile_raw}" \
#            --exclude-snps ${miss_liftover} \
#            --export vcf \
#            --out "${vcf_file}"   
# fi

# global pca plot
python "${scripts_directory}/resources/datacheck/ancestry_infer.py" \
    "${section_01_dir}/logs_b/hail.log" \
    "${vcf_file}.vcf" \
	"${genome_build}" # Use the genome build variable from the config file

# local pca plot

echo "Successfully completed script 1b"
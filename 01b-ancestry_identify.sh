#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

exec &> >(tee ${section_01c_logfile})
print_version

# 
${plink2} \
    --bfile ${bfile_raw} \
    --keep ${intersect_ids_plink} \
    --maf ${snp_maf} \
    --hwe ${snp_hwe} \
    --geno ${snp_miss} \
    --mind ${snp_imiss} \
    --make-bed \
    --out ${bfile} \
    --new-id-max-allele-len 
    --allow-extra-chr \
    --human \
    --output-chr 26 \
    --chr 1-23 \
    --threads ${nthreads}


${plink2} \
	--bfile ${bfile} \
	--set-all-var-ids @:#_\$1_\$2 \
	--make-bed \
	--out ${bfile}1 \
    --threads ${nthreads}

# calculate maf from bfile with chr:pos format
echo "Calculating MAF from bfile"
plink2 --bfile "${bfile_raw}" \
       --freq \
       --out "${section_01_dir}/"

# check against with reference panel using EasyQC
if [ -f "processed_data/genetic_data/easyQC_hrc.ecf.out" ]
then
	echo "easyqc files present from previous run which will be removed"
	rm processed_data/genetic_data/easy*
	rm processed_data/genetic_data/*easy*
else
	echo "passed file check"
fi

echo "Running global PCA"
# global pca plot
${Python_directory}python "${scripts_directory}/resources/datacheck/ancestry_infer.py" \
    "${section_01_dir}/logs_b/hail.log" \
    "${bfile_raw}" \
	"${genome_build}" \
    "${study_name}" \
    "${scripts_directory}"

# Check genome build and liftover to 38 if build is 37 using GwasDataImport
# https://github.com/MRCIEU/GwasDataImport
echo "Determining build based on reference dataset and running liftover"
${R_directory}Rscript resources/datacheck/liftover.R \
	${bfile_raw} \
	${genome_build} \
	${miss_liftover}

if [ -f ${miss_liftover} ]; then
    plink2 --bfile "${bfile_raw}" \
           --exclude-snps ${miss_liftover} \
           --make-bed \
           --out ${liftover_data}
fi


if [ "$n23" -gt "0" ]
then
    plink2 --bfile ${bfile_raw} \
	    --split-par b38 \
	    --make-bed \
	    --out ${bfile}
else
    echo "No X chromosome to split"
fi


# local pca plot to be done
echo "Successfully completed script 1b"

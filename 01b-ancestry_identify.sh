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


# adjust the easyQC script to use the correct path and file names
echo "Ancestry: ${ancestry}"

if [[ ${ancestry} == "EUR" || ${ancestry} == "AFR" || ${ancestry} == "AMR" || ${ancestry} == "EAS" || ${ancestry} == "SAS" ]]; then
    cp ${scripts_directory}/resources/genetics/1000G_${ancestry}_p3v5.TOPMed_Imputed.Allele_Freq.hg38.txt.gz ${home_directory}/processed_data/genetic_data/

else if [[ ${ancestry} == "None" ]]; then
    echo "No ancestry specified, using all population of topmed snplist and allele frequencies"
    cp ${scripts_directory}/resources/genetics/topmed.GRCh38.f8wgs.pass.nodup.mac5.maf001.tab.snplist.gz ${home_directory}/processed_data/genetic_data/
fi

replacement_text="DEFINE --pathOut "${home_directory}"/processed_data/genetic_data"
awk 'NR==3 {$0 = "'"$replacement_text"'"} 1' ${easyQCscript} > ${easyQCscript%.ecf}_edit.ecf
${R_directory}Rscript ./resources/genetics/easyQC.R ${bfile}.afreq ${easyQC} ${easyQCfile} ${easyQCscript%.ecf}_edit.ecf

rm ${home_directory}/processed_data/genetic_data/topmed.GRCh38.f8wgs.pass.nodup.mac5.maf001.tab.snplist.gz

mv ${home_directory}/processed_data/genetic_data/easyQC_hrc_edit.multi.AFCHECK.png ${home_directory}/results/02/easyQC_hrc.multi.AFCHECK.png
mv ${home_directory}/processed_data/genetic_data/easyQC_hrc_edit.rep ${home_directory}/results/02/easyQC_hrc.rep

# Remove mismatched SNPs and flip misaligned SNPs
echo "Remove mismatched SNPs and NO FLIPPING"

${plink2} \
	--bfile ${bfile} \
	--exclude ${easyQC}.mismatch_afcheck.failed.SNPs.txt \
	--make-bed \
	--out ${bfile}1 \
	--threads ${nthreads}

mv ${bfile}1.bed ${bfile}.bed
mv ${bfile}1.bim ${bfile}.bim
mv ${bfile}1.fam ${bfile}.fam




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

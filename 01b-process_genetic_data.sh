#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

exec &> >(tee ${section_01b_logfile})
print_version

# Infer genome build
echo "Inferring genome build and running liftover if necessary"
echo "Genome build in the config file is: ${genome_build}"

# Check genome build and liftover to 38 if build is 37 using GwasDataImport
# https://github.com/MRCIEU/GwasDataImport
echo "Determining build based on reference dataset and running liftover"
# Rscript will produce a map file for liftover and the SNP list that are missing from the liftover
${R_directory}Rscript resources/datacheck/liftover.R \
    ${bfile_raw} \
    ${genome_build} \
    ${miss_liftover} \
	${liftover_map} \
	${section_01_dir}

inferred_build=$(cat "${section_01_dir}/inferred_build.txt")

if [ "$inferred_build" -eq 37 ]; then
    if [ -f ${miss_liftover} ]; then
		echo "SNP missing for liftover found. Excluding them from bfile and liftovering"
        ${plink2} --bfile "${bfile_raw}" \
            --new-id-max-allele-len 70 \
            --exclude ${miss_liftover} \
            --update-map ${liftover_map} \
            --make-bed \
            --out ${bfile}
    else
		echo "No SNP missing for liftover found. Liftovering"
        ${plink2} --bfile "${bfile_raw}" \
            --new-id-max-allele-len 70 \
            --update-map ${liftover_map} \
            --make-bed \
            --out ${bfile}
    fi
elif [ "$inferred_build" -eq 38 ]; then
	# if build is 38, just copy the raw bfile to the new bfile
    ${plink2} --bfile "${bfile_raw}" \
        --new-id-max-allele-len 70 \
        --make-bed \
        --out ${bfile}
fi

# qc and format input genetic data
echo "Formatting input genetic data"
# ${intersect_ids_plink} is from 01a, it contains the IDs of samples that intersect with the genetic and methylation data

${plink2} \
    --bfile ${bfile} \
    --keep ${intersect_ids_plink} \
    --maf ${snp_maf} \
    --hwe ${snp_hwe} \
    --geno ${snp_miss} \
    --mind ${snp_imiss} \
    --make-bed \
    --out "${bfile}_format" \
    --new-id-max-allele-len 70 \
    --allow-extra-chr \
    --human \
    --output-chr 26 \
    --chr 1-23 \
    --threads ${nthreads}

# Sex check -note this is implemented in both PLINK1.9 and PLINK2

n23=`grep ^23 "${bfile}_format.bim" | wc -l`

if [ "$n23" -gt "0" ]; then
	${plink2} \
		--bfile ${bfile}_format \
		--new-id-max-allele-len 70 \
		--split-par b38 no-fail \
		--make-bed \
		--out ${bfile}_xpar_temp

	${plink2} \
		--bfile ${bfile}_xpar_temp \
		--new-id-max-allele-len 70 \
		--check-sex \
		--out ${section_01_dir}/data
	
	nprob=`grep "PROBLEM" ${section_01_dir}/data.sexcheck |wc -l`

    if [ "$nprob" -eq "0" ]
	then
		echo "There are ${nprob} individuals that failed the sex check."
	fi


	if [ "$nprob" -gt "0" ]
	then
		echo "There are ${nprob} individuals that failed the sex check."
		echo "They will be removed."
		echo "The summary is located here:"
		echo "${section_01_dir}/data.sexcheck"

		grep "PROBLEM" ${section_01_dir}/data.sexcheck | awk '{print $1, $2}' > ${bfile}_xpar_temp.failed_sexcheck
		echo "Removing individuals that failed the sex check"

		${plink2} \
			--bfile ${bfile}_xpar_temp \
			--new-id-max-allele-len 70 \
			--remove ${bfile}_xpar_temp.failed_sexcheck \
			--make-bed \
			--out ${bfile}_format

		rm ${bfile}_xpar_temp*
	fi

fi

cp ${bfile}_format.bim ${bfile}.bim.original

# format SNP ids to chr:position_A1_A2 (ascii sorted order)
echo "Formatting SNP IDs to chr:pos_A1_A2"
${plink2} \
	--bfile "${bfile}_format" \
	--new-id-max-allele-len 70 \
	--set-all-var-ids @:#_\$1_\$2 \
	--make-bed \
	--out ${bfile}1 \
    --threads ${nthreads}

cp ${bfile}_format.bim ${bfile}.bim.original1
mv ${bfile}1.bed ${bfile}.bed
mv ${bfile}1.bim ${bfile}.bim
mv ${bfile}1.fam ${bfile}.fam

touch ${SNPfail1}

${R_directory}Rscript resources/genetics/harmonization.R \
	${bfile}.bim \
	${SNPfail1}

# Checking for any duplicate SNPs
# 'exclude-mismatch': When unequal duplicate-ID variants are found, exclude every member of the group.

${plink2} \
	--bfile ${bfile} \
	--rm-dup exclude-mismatch \
	--new-id-max-allele-len 70 \
	--make-bed \
	--out ${bfile}1 \
	--threads ${nthreads}

cp ${bfile}.bim ${bfile}.bim.original2
mv ${bfile}1.bed ${bfile}.bed
mv ${bfile}1.bim ${bfile}.bim
mv ${bfile}1.fam ${bfile}.fam

# Remove SNPs with low info scores
awk '$3 < 0.80 {print $1}' <${quality_scores} > ${bfile}.lowinfoSNPs.txt

cat ${SNPfail1} ${bfile}.lowinfoSNPs.txt |sort -u >${bfile}.failed.SNPs.txt

n_failedSNPs=`wc -l ${bfile}.failed.SNPs.txt | awk '{ print $1 }'`

# Remove SNPs from data
echo "Removing ${n_failedSNPs} SNPs from data"

${plink2} \
	--bfile ${bfile} \
	--exclude ${bfile}.failed.SNPs.txt \
	--new-id-max-allele-len 70 \
	--make-bed \
	--out ${bfile}1 \
	--threads ${nthreads}

cp ${bfile}.bim ${bfile}.bim.original3
mv ${bfile}1.bed ${bfile}.bed
mv ${bfile}1.bim ${bfile}.bim
mv ${bfile}1.fam ${bfile}.fam

# Make GRMs
echo "Creating kinship matrix"
gunzip -c ${hm3_snps} > temp_hm3snps.txt

${plink2} \
	--bfile ${bfile} \
	--new-id-max-allele-len 70 \
	--extract temp_hm3snps.txt \
	--maf ${grm_maf_cutoff} \
	--make-grm-bin \
	--out ${grmfile_all} \
	--threads ${nthreads} \
	--autosome

# Create pedigree matrix if family data, otherwise remove related individuals from existing kinship and data file
if [ "${related}" = "yes" ]; then
	echo "Creating pedigree GRM"
	${R_directory}Rscript resources/relateds/grm_relateds.R ${grmfile_all} ${grmfile_relateds} ${rel_cutoff}
	
elif [ "${related}" = "no" ]; then
	echo "Removing any cryptic relateds"
	${gcta} \
		--grm ${grmfile_all} \
		--grm-cutoff ${rel_cutoff} \
		--make-grm-bin \
		--out ${grmfile_all}1
		
		mv ${grmfile_all}1.grm.N.bin ${grmfile_all}.grm.N.bin
		mv ${grmfile_all}1.grm.id ${grmfile_all}.grm.id
		mv ${grmfile_all}1.grm.bin ${grmfile_all}.grm.bin
		
		# filtering out the related samples
		${plink2} \
			--bfile ${bfile} \
			--new-id-max-allele-len 70 \
			--keep ${grmfile_all}.grm.id \
			--make-bed \
			--out ${bfile}1 \
			--threads ${nthreads}

		mv ${bfile}1.bed ${bfile}.bed
		mv ${bfile}1.bim ${bfile}.bim
		mv ${bfile}1.fam ${bfile}.fam
else 
	echo "Error: Set related flag in config to yes or no"
	exit 1
fi

#Calculate PCs
#gunzip -c ${hm3_snps_no_ld} > temp_hm3snpsnold.txt
# might to change the extract snpslist?
${plink2} \
	--bfile ${bfile} \
	--new-id-max-allele-len 70 \
	--extract temp_hm3snps.txt \
	--indep-pairwise 10000 5 0.1 \
	--maf 0.2 \
	--out ${pca} \
	--autosome \
	--threads ${nthreads}

if [ "${related}" = "no" ];then
	${plink2} \
		--bfile ${bfile} \
		--new-id-max-allele-len 70 \
		--extract ${pca}.prune.in \
		--pca 20 \
		--out ${pca} \
		--threads ${nthreads}

elif [ "${related}" = "yes" ]; then
	${plink2} \
		--bfile ${bfile} \
		--extract ${pca}.prune.in \
		--new-id-max-allele-len 70 \
		--make-bed \
		--out ${bfile}_ldpruned \
		--threads ${nthreads}

	${R_directory}Rscript resources/genetics/pcs_relateds.R \
		${bfile}_ldpruned \
		${pca} \
		${n_pcs} \
		${nthreads}
else
	echo "Error: Set related flag in config to yes or no"
	exit 1
fi

# Get genetic outliers
echo "Generating PCA plot"

${R_directory}Rscript resources/genetics/genetic_outliers.R \
	${pcs_all} \
	${pca_sd} \
	${n_pcs} \
	${genetic_outlier_ids} \
	${pcaplot}

n_outliers=`wc -l ${genetic_outlier_ids} | awk '{ print $1 }'`

if [ "${n_outliers}" -eq "0" ]
then
	echo "No genetic outliers detected"
else
	echo "There are ${n_outliers} genetic outliers detected"
	echo "They are not going to be removed from the data"
fi
# 	# Remove genetic outliers from data
# 	echo "Removing ${n_outliers} genetic outliers from data"
# 	${plink2} \
# 		--bfile ${bfile} \
# 		--new-id-max-allele-len 70 \
# 		--remove ${genetic_outlier_ids} \
# 		--make-bed \
# 		--out ${bfile}1 \
# 		--threads ${nthreads}
	
# 	mv ${bfile}1.bed ${bfile}.bed
# 	mv ${bfile}1.bim ${bfile}.bim
# 	mv ${bfile}1.fam ${bfile}.fam

# 	${gcta} \
# 		--grm ${grmfile_all} \
# 		--remove ${genetic_outlier_ids} \
# 		--make-grm-bin \
# 		--out ${grmfile_all}1 \
# 		--thread-num ${nthreads}

# mv ${grmfile_all}1.grm.N.bin ${grmfile_all}.grm.N.bin
# mv ${grmfile_all}1.grm.id ${grmfile_all}.grm.id
# mv ${grmfile_all}1.grm.bin ${grmfile_all}.grm.bin

# calculate maf from bfile with chr:pos format
echo "Calculating MAF from formatted bfile with chr:pos format"
${plink2} \
	--bfile "${bfile}" \
	--new-id-max-allele-len 70 \
	--freq \
	--out "${bfile}"

# check afreq against reference panel using EasyQC
if [ -f "processed_data/genetic_data/easyQC_topmed_edit.ecf.out" ]
then
	echo "easyqc files present from previous run which will be removed"
	rm processed_data/genetic_data/easy*
	rm processed_data/genetic_data/*easy*
else
	echo "passed easyqc file check"
fi

# Adjust the easyQC script to use the correct path and file names
echo "Ancestry: ${ancestry}"
echo "Copying reference files for EasyQC"
if [[ ${ancestry} == "EUR" || ${ancestry} == "AFR" || ${ancestry} == "AMR" || ${ancestry} == "EAS" || ${ancestry} == "SAS" ]]; then
	echo "Ancestry specified, using an imputation of ancestry-specific 1000g ref to TopMed "
    cp ${scripts_directory}/resources/genetics/1000g_${ancestry}_p3v5.topmed_imputed.maf_0.001.r2_0.3.hg38.txt.gz ${home_directory}/processed_data/genetic_data/
	replacement_text1="1000g_${ancestry}_p3v5.topmed_imputed.maf_0.001.r2_0.3.hg38.txt.gz"

elif [[ ${ancestry} == "None" ]]; then
    echo "No ancestry specified, using all population of topmed snplist and allele frequencies"
    cp ${scripts_directory}/resources/genetics/topmed.GRCh38.f8wgs.pass.nodup.mac5.maf001.tab.snplist.gz ${home_directory}/processed_data/genetic_data/
	replacement_text1=""
fi

replacement_text2="DEFINE --pathOut "${home_directory}"/processed_data/genetic_data"
awk 'NR==3 { $0 = "'"$replacement_text2"'" } 1' "${easyQCscript}" > "${easyQCscript%.ecf}_temp.ecf"

if [ -n "$replacement_text1" ]; then
	awk 'NR==30 { $0 = "\t\t--fileRef '"$replacement_text1"'" } 1' "${easyQCscript%.ecf}_temp.ecf" > "${easyQCscript%.ecf}_edit.ecf"
else
	mv "${easyQCscript%.ecf}_temp.ecf" "${easyQCscript%.ecf}_edit.ecf"
fi

rm ${easyQCscript%.ecf}_temp.ecf
easyqc_edit_ecf_cp="${genetic_processed_dir}/easyQC_topmed_edit.ecf"
mv ${easyQCscript%.ecf}_edit.ecf ${easyqc_edit_ecf_cp}

# run easyQC
echo "Running EasyQC"
${R_directory}Rscript ./resources/genetics/easyQC.R ${bfile}.afreq ${easyQC} ${easyQCfile} ${easyqc_edit_ecf_cp}

if [ -n "$replacement_text1" ]; then
	rm ${home_directory}/processed_data/genetic_data/1000g_${ancestry}_p3v5.topmed_imputed.maf_0.001.r2_0.3.hg38.txt.gz
else
	rm ${home_directory}/processed_data/genetic_data/topmed.GRCh38.f8wgs.pass.nodup.mac5.maf001.tab.snplist.gz
fi

echo "Moving allele freq check figure"

mv ${home_directory}/processed_data/genetic_data/easyQC_topmed_edit.multi.AFCHECK.png ${home_directory}/results/01/easyQC_topmed.multi.AFCHECK.png
mv ${home_directory}/processed_data/genetic_data/easyQC_topmed_edit.rep ${home_directory}/results/01/easyQC_topmed.rep

# Remove mismatched SNPs and flip misaligned SNPs
# echo "Remove mismatched SNPs and NO FLIPPING"

echo "Running global PCA"
# global pca plot using cleaned bfile
${Python_directory}python "${scripts_directory}/resources/datacheck/ancestry_infer.py" \
	"${section_01_dir}/logs_b/hail_clean.log" \
    "${bfile}" \
	"${genome_build}" \
    "${study_name}" \
    "${home_directory}" \
    "${scripts_directory}" \
	"${nthreads}" \
	"${mem}"

# From here on, we have clean data
# if [ ! "${n_outliers}" -eq "0" ]
# then

# 	echo "Recalculating PCs with outliers removed"

# 	if [ "${related}" = "no" ]
# 	then
# 		${plink2} \
# 			--bfile ${bfile} \
# 			--new-id-max-allele-len 70 \
# 			--extract ${pca}.prune.in \
# 			--pca 20 \
# 			--out ${pca} \
# 			--autosome \
# 			--threads ${nthreads}
# 	else

# 		${plink2} \
# 			--bfile ${bfile} \
# 			--new-id-max-allele-len 70 \
# 			--extract ${pca}.prune.in \
# 			--make-bed \
# 			--out ${bfile}_ldpruned \
# 			--autosome \
# 			--threads ${nthreads}

# 		${R_directory}Rscript resources/genetics/pcs_relateds.R \
# 			${bfile}_ldpruned \
# 			${pca} \
# 			${n_pcs} \
# 			${nthreads}
# 	fi

# fi

# Get frequencies, missingness, hwe, info scores
plink_files=${section_01_dir}/data*.gz
if [[ "${#plink_files[@]}" -gt 0 ]] ; then
echo "previous frequencies, missingness, hwe, info scores files present from previous run which will be removed"
	rm -f ${section_01_dir}/data.afreq.gz
	rm -f ${section_01_dir}/data.hardy.gz
	rm -f ${section_01_dir}/data.info.gz
	rm -f ${section_01_dir}/data.vmiss.gz
	rm -f ${home_directory}/processed_data/genetic_data/data.smiss.gz
else
	echo "passed file check"
fi


${plink2} \
	--bfile "${bfile}" \
	--new-id-max-allele-len 70 \
	--freq \
	--hardy \
	--missing \
	--out ${section_01_dir}/data

gzip -f -c ${quality_scores} > ${section_01_dir}/data.info.gz
gzip ${section_01_dir}/data.hardy
mv ${section_01_dir}/data.smiss ${home_directory}/processed_data/genetic_data/data.smiss
gzip ${home_directory}/processed_data/genetic_data/data.smiss
gzip ${section_01_dir}/data.vmiss
gzip ${section_01_dir}/data.afreq

# Check missingness
missingness=`zcat ${home_directory}/processed_data/genetic_data/data.smiss | awk '{ sum += $6; n++ } END { if (n > 0) print sum / n; }'`

echo "Average missingness: ${missingness}"

if (( $(bc <<< "${missingness} > 0.02") ))
then
	echo ""
	echo ""
	echo ""
	echo ""
	echo "WARNING"
	echo ""
	echo ""
	echo "Your genetic data has missingness of ${missingness}"
	echo ""
	echo "This seems high considering that you should have converted to best guess format with a very high hard call threshold"
	echo ""
	echo "Please ensure that this has been done"
fi

# Update ids
awk '{print $1,$2}' < ${bfile}.fam > ${intersect_ids_plink}
awk '{print $2}' < ${bfile}.fam > ${intersect_ids}

#rm -f ${bfile}.*~
rm temp_hm3snps.txt

# local pca plot to be done
echo "Successfully completed script 1b"
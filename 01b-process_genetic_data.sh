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
	${section_01_dir} \
	"01b"

inferred_build=$(cat "${section_01_dir}/01b_inferred_build.txt")

if [ "$inferred_build" -eq 37 ]; then
    if [ -f ${miss_liftover} ]; then
		echo "SNP missing for liftover found. Excluding them from bfile and liftovering"
        ${plink2} --bfile "${bfile_raw}" \
            --new-id-max-allele-len 70 \
            --exclude ${miss_liftover} \
            --update-map ${liftover_map} \
            --make-bed \
			--output-chr 26 \
            --out ${bfile} \
			--threads ${nthreads}
    else
		echo "No SNP missing for liftover found. Liftovering"
        ${plink2} --bfile "${bfile_raw}" \
            --new-id-max-allele-len 70 \
            --update-map ${liftover_map} \
            --make-bed \
			--output-chr 26 \
            --out ${bfile} \
			--threads ${nthreads}
    fi
elif [ "$inferred_build" -eq 38 ]; then
	# if build is 38, just copy the raw bfile to the new bfile
    ${plink2} --bfile "${bfile_raw}" \
        --new-id-max-allele-len 70 \
        --make-bed \
		--output-chr 26 \
        --out ${bfile} \
		--threads ${nthreads}
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
	${plink} \
		--bfile ${bfile}_format \
		--new-id-max-allele-len 70 \
		--split-x b38 no-fail \
		--make-bed \
		--out ${bfile}_xpar_temp \
		--threads ${nthreads}

	${plink} \
		--bfile ${bfile}_xpar_temp \
		--new-id-max-allele-len 70 \
		--check-sex \
		--out ${section_01_dir}/data \
		--threads ${nthreads}
	
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
			--output-chr 26 \
			--out ${bfile}_format \
			--threads ${nthreads}

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
	--output-chr 26 \
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
	--output-chr 26 \
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
    --sort-vars \
    --make-pgen \
    --out ${bfile_sort} \
    --output-chr 26 \
    --threads ${nthreads}

${plink2} \
    --pfile ${bfile_sort} \
    --make-bed \
    --output-chr 26 \
    --out ${bfile}1 \
    --threads ${nthreads}

rm -f ${bfile}.{pgen,psam,pvar,log}

cp ${bfile}.bim ${bfile}.bim.original3
mv ${bfile}1.bed ${bfile}.bed
mv ${bfile}1.bim ${bfile}.bim
mv ${bfile}1.fam ${bfile}.fam

# =====================================================================
# Make GRMs  (two parallel tracks: GCTA + KING/PC-Relate)
# =====================================================================

gunzip -c ${hm3_snps} > temp_hm3snps.txt

# ====================================================================
# TRACK 1: GCTA dense kinship matrix (on full sample, pre-filter)
# ====================================================================
echo "Creating kinship matrix (GCTA)"

${plink2} \
    --bfile ${bfile} \
    --new-id-max-allele-len 70 \
    --extract temp_hm3snps.txt \
    --maf ${grm_maf_cutoff} \
    --make-grm-bin \
    --out ${grmfile_all} \
    --threads ${nthreads} \
    --autosome

echo "Sample size in creating GCTA kinship matrix: $(wc -l < "${grmfile_all}.grm.id")"

# # =====================================================================
# # Make GRMs and plot distributions
# # =====================================================================
# # Pre-filter kinship calculation for plotting/QC purposes
echo "Preparing KING input bfile"

${plink} \
    --bfile ${bfile} \
    --new-id-max-allele-len 70 \
    --make-bed \
    --exclude range "${exclude_highld_region}" \
    --out "${bfile}_king_input" \
    --output-chr 26 \
    --autosome \
    --threads ${nthreads}

${king} \
    -b ${bfile}_king_input.bed \
    --kinship \
    --cpus ${nthreads} \
    --prefix ${grmfile_king}_prefilter

echo "Plotting prefilter-GRM distributions (GCTA + KING) for full sample"
${R_directory}Rscript resources/relateds/gcta_king_grm_distri.R \
    "${grmfile_all}" \
    "${rel_cutoff}" \
    "${grm_distribution}_01b" \
    "${grmfile_king}_prefilter.kin0" \
    "king"

# Derive KING degree from rel_cutoff (GCTA scale: kinship = rel_cutoff / 2)
kinship_cutoff=$(echo "${rel_cutoff} / 2" | bc -l)
if   (( $(echo "$kinship_cutoff > 0.0884" | bc -l) )); then king_degree=2
elif (( $(echo "$kinship_cutoff > 0.0442" | bc -l) )); then king_degree=3
else king_degree=4
fi
echo ">> KING degree derived from rel_cutoff (${rel_cutoff}): --degree ${king_degree}"

# =====================================================================
# PRIMARY SPLIT: structured flag from config determines kinship method
# =====================================================================

if [ "${structured}" = "yes" ]; then
    echo ">> Structured population: using PC-Relate for kinship estimation"

    if [ "${related}" = "yes" ]; then
        echo ">> SCENARIO 1: Related + Structured -> PC-Relate sparse GRM (keep relatives)"
        pcrelate_mode="keep"
    elif [ "${related}" = "no" ]; then
        echo ">> SCENARIO 3: Unrelated + Structured -> PC-Relate remove cryptic relatedness"
        pcrelate_mode="remove"
    else
        echo "Error: Please set related flag to 'yes' or 'no' in config."
        exit 1
    fi

    ${R_directory}Rscript resources/relateds/compute_pcrelate_grm.R \
        "${bfile}_king_input" \
        "${grmfile_pcrelate}" \
        "${n_pcs}" \
        "${nthreads}" \
        "${rel_cutoff}" \
        "${pcrelate_mode}"

    ${R_directory}Rscript resources/relateds/gcta_king_grm_distri.R \
        "${grmfile_all}" \
        "${rel_cutoff}" \
        "${grm_distribution}_01b_pcrelate" \
        "${grmfile_pcrelate}.pcrelate.kin0.txt" \
        "pc_relate"

    if [ "${related}" = "no" ]; then
        n_before=$(wc -l < "${bfile}.fam")
        n_remove=$(wc -l < "${grmfile_pcrelate}.remove_ids.txt")
        echo ">> Sample size before PC-Relate removal: ${n_before}"
        echo ">> Samples to remove (cryptic relateds): ${n_remove}"

        ${plink2} \
            --bfile ${bfile} \
            --remove ${grmfile_pcrelate}.remove_ids.txt \
            --new-id-max-allele-len 70 \
            --output-chr 26 \
            --make-bed --out ${bfile}_unrel_final --threads ${nthreads}

        mv ${bfile}_unrel_final.bed ${bfile}.bed
        mv ${bfile}_unrel_final.bim ${bfile}.bim
        mv ${bfile}_unrel_final.fam ${bfile}.fam

        n_after=$(wc -l < "${bfile}.fam")
        echo ">> Sample size after PC-Relate removal: ${n_after}"
    fi

    # PCs already computed by compute_pcrelate_grm.R (PC-Air)
    cp "${grmfile_pcrelate}.pca.eigenvec" "${pca}.eigenvec"

elif [ "${structured}" = "no" ]; then
    echo ">> Non-structured population: using KING for kinship estimation"

    if [ "${related}" = "yes" ]; then
        echo ">> SCENARIO 2: Related + Non-Structured -> KING sparse GRM"

        ${king} \
            -b ${bfile}_king_input.bed \
            --related \
            --degree ${king_degree} \
            --cpus ${nthreads} \
            --prefix ${grmfile_king}_filter

        if [ -f "${grmfile_king}_filter.kin0" ]; then
            awk '{ print $2, $4, $14 }' ${grmfile_king}_filter.kin0 | grep -v "4th" | sed 1d > ${grmfile_king}_filter.kin0.formatted

            ${R_directory}Rscript resources/relateds/pedFAM.R \
                ${bfile}_king_input.fam \
                ${grmfile_king}_filter.kin0.formatted \
                ${grmfile_king}_sparse
        fi

    elif [ "${related}" = "no" ]; then
        echo ">> SCENARIO 4: Unrelated + Non-Structured -> KING remove cryptic relatedness"
        echo ">> rel_cutoff: ${rel_cutoff} -> KING --degree ${king_degree}"

        ${king} \
            -b ${bfile}_king_input.bed \
            --unrelated \
            --degree ${king_degree} \
            --cpus ${nthreads} \
            --prefix ${grmfile_king}_filter

        n_before=$(wc -l < "${bfile}.fam")
        echo ">> Sample size before KING removal: ${n_before}"

        echo "Filtering genotypes to keep only the unrelated list"
        ${plink2} \
            --bfile ${bfile} \
            --keep ${grmfile_king}_filterunrelated.txt \
            --new-id-max-allele-len 70 \
            --output-chr 26 \
            --make-bed --out ${bfile}_unrel_final --threads ${nthreads}
        mv ${bfile}_unrel_final.bed ${bfile}.bed
        mv ${bfile}_unrel_final.bim ${bfile}.bim
        mv ${bfile}_unrel_final.fam ${bfile}.fam

        n_after=$(wc -l < "${bfile}.fam")
        echo ">> Sample size after KING removal: ${n_after}"

    else
        echo "Error: Please set related flag to 'yes' or 'no' in config."
        exit 1
    fi

    # PC calculation (structured=no only; structured=yes, then PCs already done by PC-Air)
    ${plink2} \
        --bfile ${bfile} \
        --new-id-max-allele-len 70 \
        --extract temp_hm3snps.txt \
        --indep-pairwise 10000 5 0.1 \
        --maf 0.2 \
        --out ${pca} \
        --autosome \
        --threads ${nthreads}

    if [ "${related}" = "no" ]; then
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
    fi

else
    echo "Error: Please set structured flag to 'yes' or 'no' in config."
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
	echo "They are not going to be removed from the sample in QC stage 1"
fi

# calculate maf from bfile with chr:pos format
echo "Calculating MAF from formatted bfile with chr:pos format"
${plink2} \
	--bfile "${bfile}" \
	--new-id-max-allele-len 70 \
	--output-chr 26 \
	--freq \
	--out "${bfile}" \
	--threads ${nthreads}

# check afreq against reference panel using EasyQC
if [ -f "${home_directory}/processed_data/genetic_data/easyQC_topmed_edit.ecf.out" ]
then
	echo "easyqc files present from previous run which will be removed"
	rm ${home_directory}/processed_data/genetic_data/easy*
	rm ${home_directory}/processed_data/genetic_data/*easy*
else
	echo "passed easyqc file check"
fi

# Replace indels with I&D for easyQC; other failed SNPs have been removed already in line 145
echo "Replacing indels with I&D for easyQC"

${R_directory}Rscript resources/genetics/harmonization_indel.R \
	"${bfile}.afreq" \
	"${bfile}.easyQC.afreq"

# Adjust the easyQC script to use the correct path and file names
echo "Ancestry: ${ancestry}"
echo "Copying reference files for EasyQC"
if [[ ${ancestry} == "EUR" || ${ancestry} == "AFR" || ${ancestry} == "AMR" || ${ancestry} == "EAS" || ${ancestry} == "SAS" ]]; then
	echo "Ancestry specified, using an imputation of ancestry-specific 1000g ref to TopMed "
    cp ${scripts_directory}/resources/genetics/1000g_${ancestry}_p3v5.topmed_imputed.maf_0.001.r2_0.3.indelrecoded.hg38.txt.gz ${home_directory}/processed_data/genetic_data/
	replacement_text1="1000g_${ancestry}_p3v5.topmed_imputed.maf_0.001.r2_0.3.indelrecoded.hg38.txt.gz"

elif [[ ${ancestry} == "None" ]]; then
    echo "No ancestry specified, using all population of topmed snplist and allele frequencies"
    cp ${scripts_directory}/resources/genetics/topmed.GRCh38.f8wgs.pass.nodup.mac5.maf001.indelrecoded.tab.snplist.gz ${home_directory}/processed_data/genetic_data/
	replacement_text1=""
fi

replacement_text2="DEFINE --pathOut "${home_directory}"/processed_data/genetic_data"
awk 'NR==3 { $0 = "'"$replacement_text2"'" } 1' "${easyQCscript}" > "${easyQCscript%.ecf}_temp.ecf"

if [ -n "$replacement_text1" ]; then
	awk 'NR==30 { $0 = "\t\t--fileRef '"$replacement_text1"'" } 1' "${easyQCscript%.ecf}_temp.ecf" > "${easyQCscript%.ecf}_edit.ecf"
else
	mv "${easyQCscript%.ecf}_temp.ecf" "${easyQCscript%.ecf}_edit.ecf"
fi

rm -f ${easyQCscript%.ecf}_temp.ecf
easyqc_edit_ecf_cp="${genetic_processed_dir}/easyQC_topmed_edit.ecf"
mv ${easyQCscript%.ecf}_edit.ecf ${easyqc_edit_ecf_cp}

# run easyQC
echo "Running EasyQC"
${R_directory}Rscript ./resources/genetics/easyQC.R ${bfile}.easyQC.afreq ${easyQC} ${easyQCfile} ${easyqc_edit_ecf_cp}

if [ -n "$replacement_text1" ]; then
	rm ${home_directory}/processed_data/genetic_data/1000g_${ancestry}_p3v5.topmed_imputed.maf_0.001.r2_0.3.indelrecoded.hg38.txt.gz
else
	rm ${home_directory}/processed_data/genetic_data/topmed.GRCh38.f8wgs.pass.nodup.mac5.maf001.indelrecoded.tab.snplist.gz
fi

echo "Moving allele freq check figure"

mv ${home_directory}/processed_data/genetic_data/easyQC_topmed_edit.multi.AFCHECK.png ${home_directory}/results/01/easyQC_topmed.multi.AFCHECK.png
mv ${home_directory}/processed_data/genetic_data/easyQC_topmed_edit.rep ${home_directory}/results/01/easyQC_topmed.rep
cp ${home_directory}/processed_data/genetic_data/data.easyqc.AFCHECK.outlier.txt ${home_directory}/results/01/data.easyqc.AFCHECK.outlier.txt

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
	--output-chr 26 \
	--out ${section_01_dir}/data \
	--threads ${nthreads}

gzip -f -c ${quality_scores} > ${section_01_dir}/data.info.gz
gzip -f ${section_01_dir}/data.hardy
gzip -f ${section_01_dir}/data.smiss
gzip -f ${section_01_dir}/data.vmiss
gzip -f ${section_01_dir}/data.afreq

echo "Moving smiss file to processed_data/genetic_data"
mv ${section_01_dir}/data.smiss.gz ${home_directory}/processed_data/genetic_data/data.smiss.gz

# Check missingness
# zcat auto-resolve data.smiss -> data.smiss.gz on Linux; keep ".gz" explicit for readability.
missingness=`zcat ${home_directory}/processed_data/genetic_data/data.smiss.gz | awk '{ sum += $6; n++ } END { if (n > 0) print sum / n; }'`

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
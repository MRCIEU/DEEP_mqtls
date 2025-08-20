#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

exec &> >(tee ${section_01c_logfile})
print_version

# Check phenotype distributions, remove outliers
echo "Check phenotypes"
# smoking from cohorts
${R_directory}Rscript resources/datacheck/check_phenotypes.R \
		${betas} \
		${methylation_no_outliers} \
		${cohort_descriptives_commonids} \
		${methylation_summary} \
		${intersect_ids} \
		${covariates} \
		${covariates_intersect} \
		${bfile}.fam \
		${bfile}.bim

# Predict age and smoking
echo "predict age and smoking"
${R_directory}Rscript resources/datacheck/predict_age_smoking.R \
		${methylation_no_outliers} \
		${bfile}.fam \
		${smoking_pred} \
		${smoking_pred_plot} \
		${smoking_pred_SD} \
		${covariates} 

# cell counts and correlations
echo "cell counts and correlations"

${R_directory}Rscript resources/cellcounts/cellcounts_epiDISH.R \
		${tissue} \
		${age} \
		${methylation_array} \
        ${methylation_no_outliers} \
        ${cellcounts_cov} \
        ${cellcounts_plot} \
        ${cellcounts_summary} \
		${scripts_directory} 

if [ "${measured_cellcounts}" != "NULL" ] && [ -f "${measured_cellcounts}" ]; then
	echo "Comparing measured with predicted cellcounts"
    ${R_directory}Rscript resources/cellcounts/correlation.R \
		${cellcounts_cov} \
		${measured_cellcounts} \
		${cor_matrix} \
		${cor_plot} \
		${scripts_directory}
elif [ "${measured_cellcounts}" == "NULL" ]; then
	echo "No measured cell counts available for comparison; only compare predicted cell counts if multiple reference available"
	${R_directory}Rscript resources/cellcounts/correlation.R \
		${cellcounts_cov} \
		0 \
		${cor_matrix} \
		${cor_plot} \
		${scripts_directory}
else
	echo "Error: measured_cellcounts is not 'NULL' but the file does not exist. Please check your config file."
	exit 1
fi

# EWAS of age and smoking
echo "EWAS of age and smoking"
${R_directory}Rscript resources/datacheck/ewas_age_smoking.R \
        ${methylation_no_outliers} \
        ${cellcounts_cov} \
        ${cellcounts_plot} \
        ${cellcounts_summary}

# Combine covariates
# smoking from prediction
echo "Combining covariates for mQTL analysis"
${R_directory}Rscript resources/genetics/covariates.R \
	${covariates_intersect} \
	${pcs_all} \
	${cellcounts_cov} \
	${smoking_pred}.txt \
	${bfile}.fam \
	${covariates_combined}

echo "Successfully completed script 1c"
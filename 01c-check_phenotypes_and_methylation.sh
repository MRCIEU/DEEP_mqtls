#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

exec &> >(tee ${section_01c_logfile})
print_version

containsElement () {
	local e
	for e in "${@:2}"; do [[ "$e" == "$1" ]] && return 0; done
	echo "There is no method for ${1}."
	echo "Please run:"
	echo "./01c-check_phenotypes_and_methylation.sh [arg]"
	echo "where arg is an optional argument that can be one of:"
	printf '%s\n' ${@:2}
	return 1
}

arg="all"
declare -a sections=('all' 'methy_outlier' 'check_phenotype' 'predict_age_smoking' 'cell_counts' 'ewas' 'combine_covariates')

if [ -n "${1}" ]; then
	arg="${1}"
	containsElement ${1} ${sections[@]}
fi

section_message () {

	echo "-----------------------------------------------"
	echo ""
	echo "$1 section"
	echo ""
	echo "to run this part on its own type:"
	echo "$ ./01c-check_phenotypes_and_methylation.sh $1"
	echo ""
	echo "-----------------------------------------------"
	echo ""
	echo ""

}

if [ "$arg" = "methy_outlier" ] || [ "$arg" = "all" ]
then
	section_message "methy_outlier"

	echo "Removing methylation outliers"
	${R_directory}Rscript resources/methylation/remove_outliers.R \
		${betas} \
		${methylation_no_outliers} \
		${cohort_descriptives_commonids} \
		${methylation_summary} \
		${intersect_ids} \
		${covariates} \
		${covariates_intersect} \
		${bfile}.fam \
		${bfile}.bim

fi

if [ "$arg" = "check_phenotype" ] || [ "$arg" = "all" ]
then
	section_message "check_phenotype"

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

fi

if [ "$arg" = "predict_age_smoking" ] || [ "$arg" = "all" ]
then
	section_message "predict_age_smoking"

	echo "Predict age and smoking"
	${R_directory}Rscript resources/datacheck/predict_age_smoking.R \
		${methylation_no_outliers} \
		${bfile}.fam \
		${smoking_pred} \
		${smoking_pred_plot} \
		${smoking_pred_SD} \
		${covariates}

fi

if [ "$arg" = "cell_counts" ] || [ "$arg" = "all" ]
then
	section_message "cell_counts"

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

fi

if [ "$arg" = "ewas" ] || [ "$arg" = "all" ]
then
	section_message "ewas"

	echo "EWAS of age and smoking"
	${R_directory}Rscript resources/datacheck/ewas_age_smoking.R \
		${methylation_no_outliers} \
		${cellcounts_cov} \
		${cellcounts_plot} \
		${cellcounts_summary}

fi

if [ "$arg" = "combine_covariates" ] || [ "$arg" = "all" ]
then
	section_message "combine_covariates"

	echo "Combining covariates for mQTL analysis"
	${R_directory}Rscript resources/genetics/covariates.R \
		${covariates_intersect} \
		${pcs_all} \
		${cellcounts_cov} \
		${smoking_pred}.txt \
		${bfile}.fam \
		${covariates_combined}

fi

if [ "$arg" = "all" ]
then
	echo "Successfully completed script 01c"
fi
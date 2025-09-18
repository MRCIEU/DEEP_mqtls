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
declare -a sections=('all' 'methy_outlier' 'check_phenotype' 'predict_age_smoking' 'cell_counts' 'ewas' 'meth_pcs' 'combine_covariates')

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
		${intersect_ids} \
		${meth_ids} \
		${covariates} \
		${bfile}.fam \
		${bfile}.bim \
		${methylation_no_outliers_gwas} \
		${methylation_no_outliers_ewas} \
		${cohort_descriptives_was} \
		${methylation_summary_gwas} \
		${methylation_summary_ewas} \
		${covariates_intersect_gwas} \
		${covariates_intersect_ewas} \

fi

if [ "$arg" = "check_phenotype" ] || [ "$arg" = "all" ]
then
	section_message "check_phenotype"

	echo "Check phenotypes"
	# smoking from cohorts
	${R_directory}Rscript resources/datacheck/check_phenotypes.R \
		${phenotypes} \
		${bfile}.fam \
		${meth_ids} \
		${raw_phenotype_distribution_plot} \
		${raw_phenotype_summary_file} \
		${study_name} \
		${edited_phenotype_distribution_plot} \
		${edited_phenotype_summary_file} \
		${winsorized_phenotype_file} 

fi

if [ "$arg" = "predict_age_smoking" ] || [ "$arg" = "all" ]
then
	section_message "predict_age_smoking"

	echo "Predict age and smoking"
	${R_directory}Rscript resources/datacheck/predict_age_smoking.R \
		${methylation_no_outliers_ewas} \
		${winsorized_phenotype_file} \
		${bfile}.fam \
		${study_name} \
		${age_smoking_prediction_plot} \
		${updated_phenotype_file} \
		${predicted_smoking}

fi

if [ "$arg" = "cell_counts" ] || [ "$arg" = "all" ]
then
	section_message "cell_counts"

	
	echo "cell counts and correlations"

	${R_directory}Rscript resources/cellcounts/cellcounts_epiDISH.R \
		"${tissue}" \
		"${age}" \
		"${methylation_array}" \
		"${methylation_no_outliers_gwas}" \
		"${cellcounts_cov}" \
		"${cellcounts_plot}" \
		"${cellcounts_summary}" \
		"${scripts_directory}"

	if [ "${measured_cellcounts}" != "NULL" ] && [ -f "${measured_cellcounts}" ]; then
		echo "Comparing measured cellcounts with predicted ones"
		
		${R_directory}Rscript resources/cellcounts/correlation.R \
			"${cellcounts_cov}" \
			"${format_measured_cellcounts}" \
			"${cellcounts_cov_total}" \
			"${cor_matrix}" \
			"${cor_plot_ori}" \
			"${cor_plot_comb}" \
			"${study_name}" \
			"${home_directory}"

	elif [ "${measured_cellcounts}" == "NULL" ]; then
		echo "No measured cell counts available for comparison; only compare predicted cell counts if multiple reference available"
		${R_directory}Rscript resources/cellcounts/correlation.R \
			"${cellcounts_cov}" \
			0 \
			"${cellcounts_cov_total}" \
			"${cor_matrix}" \
			"${cor_plot_ori}" \
			"${cor_plot_comb}" \
			"${study_name}" \
			"${home_directory}"
	else
		echo "Error: measured_cellcounts is not 'NULL' but the file does not exist. Please check your config file."
		exit 1
	fi

fi

if [ "$arg" = "ewas" ] || [ "$arg" = "all" ]
then
	section_message "ewas"

	echo "EWAS of age and smoking - cell count adjusted"
	${R_directory}Rscript resources/datacheck/ewas_age_smoking.R \
		${methylation_no_outliers_ewas} \
		${updated_phenotype_file} \
		${methylation_array} \
		${cellcounts_cov} \
		${study_name} \
		${qc1_ewas_stats} \
		${qc1_ewas_report} \
		${scripts_directory} \
		${study_specific_vars}

fi

if [ "$arg" = "meth_pcs" ] || [ "$arg" = "all" ]
then
	section_message "meth_pcs"

	echo "meth PCs"
	${R_directory}Rscript resources/datacheck/methylation_pc_check.R \
		${methylation_no_outliers_ewas} \
		${updated_phenotype_file} \
		${cellcounts_cov} \
		${study_name} \
		${pcs_all} \
		${meth_pcs_scree_plot} \
		${meth_pcs_screePC1PC2_plot} \
		${meth_pcs_screePC3PC4_plot} \
		${pc_var_association_plot} \
		${scripts_directory} \
		${study_specific_vars}

fi

if [ "$arg" = "combine_covariates" ] || [ "$arg" = "all" ]
then
	section_message "combine_covariates"

	echo "Combining covariates for mQTL analysis"
	${R_directory}Rscript resources/genetics/covariates.R \
		"${covariates_intersect_gwas}" \
		"${pcs_all}" \
		"${cellcounts_cov}" \
		"${predicted_smoking}" \
		"${bfile}.fam" \
		"${covariates_combined}"

fi

if [ "$arg" = "all" ]
then
	echo "Successfully completed script 01c"
fi

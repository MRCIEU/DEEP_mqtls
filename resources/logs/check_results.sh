#!/usr/bin/env bash

check_results_01a () {

	if [ -f "${cohort_descriptives}" ]; then
		echo "Cohort descriptives file present"
	else
		echo "Cohort descriptives file is absent."
		exit 1
	fi

	if [ -f "${ageplot}" ]; then
		echo "Age distribution plot file present"
	else
		echo "Age distribution plot file is absent."
		exit 1
	fi

	if [ -f "${snpchrtxt}" ] && [ -f "${snpchrplot}" ]; then
		echo "SNP chromosome distribution files present"
	else
		echo "SNP chromosome distribution files are absent."
		exit 1
	fi

	if [ -f "${quality_scores_plot}" ]; then
		echo "Quality scores plot file present"
	else
		echo "Quality scores plot file is absent."
		exit 1
	fi

	if [ -f "${sex_pred_plot}" ]; then
		echo "Sex prediction plot file present"
	else
		echo "The sex prediction plot file is missing. Please check if the samples are of a single sex (e.g., all male or all female) or if the sex probes are missing from the dataset."
	fi
}

check_results_01b () {

	if [ -f "${section_01_dir}/01b_inferred_build.txt" ]; then
		echo "Inferred build file from 01b present"
	else
		echo "Problem: inferred build file from 01b is absent"
		exit 1
	fi

	if [ -f "${section_01_dir}/pcaplot.pdf" ]; then
		echo "PCA plot present"
	else
		echo "Problem: PCA plot is absent"
		exit 1
	fi

	if [ -f "${section_01_dir}/${study_name}_globalPCA.png" ]; then
		echo "Global PCA plot present"
	else
		echo "Problem: Global PCA plot is absent"
		exit 1
	fi

	if [ -f "${section_01_dir}/easyQC_topmed.multi.AFCHECK.png" ] && [ -f "${section_01_dir}/easyQC_topmed.rep" ]; then
    	echo "easyQC plot and results present"
	else
    	echo "Problem: easyQC plot or results are absent"
    	exit 1
	fi

	if [ -f "${section_01_dir}/data.afreq.gz" ]; then
		echo "Allele frequency file present"
	else
		echo "Problem: Allele frequency file is absent"
		exit 1
	fi

	if [ -f "${section_01_dir}/data.hardy.gz" ]; then
		echo "HWE file present"
	else
		echo "Problem: HWE file is absent"
		exit 1
	fi

	if [ -f "${section_01_dir}/data.info.gz" ]; then
		echo "Imputation quality file present"
	else
		echo "Problem: Imputation quality file is absent"
		exit 1
	fi

	if [ -f "${home_directory}/processed_data/genetic_data/data.smiss.gz" ]; then
		echo "Sample missingness file present in processed_data folder"
	else
		echo "Problem: Sample missingness file is absent"
		exit 1
	fi
	
	# confirm that sample missingness file is not in section_01_dir
	if [ -f "${section_01_dir}/data.smiss.gz" ]; then
    	echo "ERROR: Sample missingness file 'data.smiss.gz' is still in ${section_01_dir}."
    	echo "Please move it to 'processed_data/genetic_data' and delete it from the results folder manually."
    	exit 1
	fi

}

check_results_01c () {

	if [ -f "${section_01_dir}/methylation_summary.RData" ] && [ -f "${section_01_dir}/methylation_summary_gwas.RData" ] && [ -f "${section_01_dir}/methylation_summary_ewas.RData" ]; then
		echo "Three methylation_summary RData files are present"
	else
		echo "Problem: one or more methylation_summary RData files are absent"
		exit 1
	fi

    if [ -f "${section_01_dir}/cohort_descriptives_was.RData" ]; then
		echo "cohort_descriptives_was.RData is present"
	else
		echo "Problem: cohort_descriptives_was.RData is absent"
		exit 1
	fi

	if [ -f "${section_01_dir}/cellcounts_summary.txt" ]; then
		echo "Summary statistics of cell counts are present"
	else
		echo "Problem: summary statistics of cell counts are absent"
		exit 1
	fi

	if [ -f "${section_01_dir}/cellcounts_plot.pdf" ]; then
		echo "Plots of cell counts are present"
	else
		echo "Problem: plots of cell counts are absent"
		exit 1
	fi

	if [ ! -f "${section_01_dir}/cor_plot_ori.pdf" ]; then
    	echo "Problem: ori correlation plot is absent"
    	exit 1
	fi
	if [ ! -f "${section_01_dir}/cor_plot_comb.pdf" ]; then
    	echo "Problem: combined correlation plot is absent"
    	exit 1
	fi
	echo "Ori and combined correlation plots of cell counts are present"

	if [  -f "${section_01_dir}/cor_matrix.txt" ]; then
		echo "Correlation matrix of observed vs predicted cell counts is present"
	else
		echo "Problem: correlation matrix of observed vs predicted cell counts is absent"
		exit 1
	fi

	if [ -d "${section_01_dir}/cellcounts_comp" ]; then
    pdf_count=$(find "${section_01_dir}/cellcounts_comp" -name "*.pdf" | wc -l)
    	if [ $pdf_count -gt 0 ]; then
        	echo "cellcounts_comp folder contains $pdf_count PDF files"
    	else
        	echo "cellcounts_comp folder contains no PDF files"
   		fi
	fi

	if [ -f "${age_smoking_prediction_plot}"*.jpg ]; then
		echo "Smoking and age prediction plot is present"
	else
		echo "Problem: Smoking and age prediction plot file not present"
		exit 1
	fi

	if ls "${meth_pcs_scree_plot}"*.pdf >/dev/null 2>&1 && ls "${meth_pcs_PC1PC2_plot}"*.pdf >/dev/null 2>&1 && ls "${meth_pcs_PC3PC4_plot}"*.pdf >/dev/null 2>&1; then
    	echo "The meth PCs plots are present"
	else
    	echo "Problem: One or more of the meth PCs plots are absent"
    	exit 1
	fi

	if ls "${section_01_dir}/pc_var_association_plot_${study_name}"*.pdf >/dev/null 2>&1; then
    	echo "The pc_var_association plot is present"
	else
    	echo "Problem: The pc_var_association plot is absent"
    	exit 1
	fi

	shopt -s nullglob
    ewas_stats_files=( "${qc1_ewas_stats}"* )
    shopt -u nullglob
    if [ ${#ewas_stats_files[@]} -gt 0 ]; then
        echo "EWAS stats files present (${#ewas_stats_files[@]}):"
        for f in "${ewas_stats_files[@]}"; do
            printf '  %s\n' "${f##*/}"
        done
    else
        echo "Problem: EWAS stats file is absent"
        exit 1
    fi

    shopt -s nullglob
    ewas_report_files=( "${qc1_ewas_report}"* )
    shopt -u nullglob
    if [ ${#ewas_report_files[@]} -gt 0 ]; then
        echo "EWAS report files present (${#ewas_report_files[@]}):"
        for f in "${ewas_report_files[@]}"; do
            printf '  %s\n' "${f##*/}"
        done
    else
        echo "Problem: EWAS report file is absent"
        exit 1
    fi

	shopt -s nullglob
    cell_count_scree_files=( "${cell_count_scree_plot}"* )
    shopt -u nullglob
    if [ ${#cell_count_scree_files[@]} -gt 0 ]; then
        echo "Cell count scree files present (${#cell_count_scree_files[@]}):"
        for f in "${cell_count_scree_files[@]}"; do
            printf '  %s\n' "${f##*/}"
        done
    else
        echo "Problem: Cell count scree file is absent"
        exit 1
    fi

	if [ -f "${section_01_dir}/raw_phenotype_distribution_plot_${study_name}.jpg" ] && [ -f "${section_01_dir}/raw_phenotypes_summary_${study_name}.Rdata" ]; then
    	echo "Raw phenotypes summary files are present"
	else
    	echo "Problem: One or more raw phenotypes summary files are absent"
    	exit 1
	fi
	
	# Check if winsorisation was applied to numeric phenotypes in input data
	if grep -q "no numeric vars to winsorise" "${section_01_dir}/logs_c"/log*.txt; then
    	echo "Winsorisation was not applied to phenotype input data so no edited phenotypes output files expected"
	else
    # Winsorisation was applied to numeric phenotypes in input data
    	if [ -f "${section_01_dir}/edited_phenotypes_summary_${study_name}.Rdata" ] && [ -f "${section_01_dir}/edited_phenotype_distribution_plot_${study_name}.jpg" ]; then
        	echo "Edited phenotypes summary files are present"
    	else
        echo "Problem: One or more edited phenotypes summary files are absent"
        exit 1
  	  fi
fi
}

check_results_01d () {

    if [ -f "${methylation_processed_dir}/mqtl_pos_ctr_1d_filt.tsv" ]; then
        echo "Filtered positive control file is present in processed_data folder"
    else
        echo "Problem: filtered positive control file is absent in processed_data folder"
        exit 1
    fi

    if [ ! -d "${section_01_dir}/01d" ]; then
        echo "Problem: directory ${section_01_dir}/01d does not exist"
        exit 1
    fi

    local suffixes=(
        ".fastGWA.gz"
        "_manhattan.pdf"
        "_nocisChr_manhattan.pdf"
        "_qqplot.jpeg"
        "_nocisChr_qqplot.jpeg"
    )

    while IFS=$'\t' read -r positive_control_cpg _; do
        [ -z "${positive_control_cpg}" ] && continue

        echo "Checking positive control: ${positive_control_cpg}"

        ## ---- Positive control ----
        pos_base="${section_01_dir}/01d/positive_control_untransformed_${positive_control_cpg}"

        for suffix in "${suffixes[@]}"; do
            file="${pos_base}${suffix}"

            if [ ! -f "$file" ]; then
                echo "Problem: positive control ${positive_control_cpg} missing file $file"
                exit 1
            fi

            case "$suffix" in
                ".fastGWA.gz")
                    echo "positive control ${positive_control_cpg} results present"
                    ;;
                "_manhattan.pdf")
                    echo "positive control ${positive_control_cpg} Manhattan plot present"
                    ;;
                "_nocisChr_manhattan.pdf")
                    echo "positive control ${positive_control_cpg} no cis chromosome Manhattan plot present"
                    ;;
                "_qqplot.jpeg")
                    echo "positive control ${positive_control_cpg} QQ plot present"
                    ;;
                "_nocisChr_qqplot.jpeg")
                    echo "positive control ${positive_control_cpg} no cis chromosome QQ plot present"
                    ;;
            esac
        done

        ## ---- Negative control ----
        neg="NEG_${positive_control_cpg}"
        echo "Checking negative control: ${neg}"

        neg_base="${section_01_dir}/01d/negative_control_untransformed_${neg}"

        for suffix in "${suffixes[@]}"; do
            file="${neg_base}${suffix}"

            if [ ! -f "$file" ]; then
                echo "Problem: negative control ${neg} missing file $file"
                exit 1
            fi

            case "$suffix" in
                ".fastGWA.gz")
                    echo "negative control ${neg} results present"
                    ;;
                "_manhattan.pdf")
                    echo "negative control ${neg} Manhattan plot present"
                    ;;
                "_nocisChr_manhattan.pdf")
                    echo "negative control ${neg} no cis chromosome Manhattan plot present"
                    ;;
                "_qqplot.jpeg")
                    echo "negative control ${neg} QQ plot present"
                    ;;
                "_nocisChr_qqplot.jpeg")
                    echo "negative control ${neg} no cis chromosome QQ plot present"
                    ;;
            esac
        done

    done < <(
        awk -F'\t' '
            NR > 1 {
                gsub(/^[ \t]+|[ \t]+$/, "", $1)
                if ($1 != "" && $1 !~ /^#/)
                    print $1
            }
        ' "${methylation_processed_dir}/mqtl_pos_ctr_1d_filt.tsv"
    )

}



check_results_01e () {

    if [ ! -d "${section_01_dir}/01e" ]; then
        echo "Problem: directory ${section_01_dir}/01e does not exist"
        exit 1
    fi

    for i in {1..5}; do
        pc="PC${i}"

        if [ -f "${section_01_dir}/01e/gwas_${pc}.fastGWA.gz" ]; then
            echo "GWAS ${pc} results present"
        else
            echo "Problem: GWAS ${pc} results file not present: ${section_01_dir}/01e/gwas_${pc}.fastGWA.gz"
            exit 1
        fi

        if [ -f "${section_01_dir}/01e/gwas_${pc}_manhattan_beta.pdf" ]; then
            echo "GWAS ${pc} beta Manhattan plot present"
        else
            echo "Problem: GWAS ${pc} beta Manhattan plot file not present: ${section_01_dir}/01e/gwas_${pc}_manhattan_beta.pdf"
            exit 1
        fi
    done
}

check_results_01f () {

	if [ -f "${section_01_dir}/${study_name}.qc_objects_shrink.rda" ]; then
		echo "qc_objects_shrink R object present"
	else
		echo "WARNING: qc_objects_shrink R object is absent. We won't be able to perform cross-cohort normalisation for you. Please confirm if you perform normalization by meffil or not."
	fi

}

check_results_01g () {

    if [ ! -f "${methylation_processed_dir}/mqtl_pos_ctr_1g_filt.tsv" ]; then
        echo "Problem: filtered positive control file is absent in processed_data folder"
        exit 1
    fi

    if [ ! -d "${section_01_dir}/01g" ]; then
        echo "Problem: directory ${section_01_dir}/01g does not exist"
        exit 1
    fi

    if [ -f "${section_01_dir}/01g_inferred_build.txt" ]; then
		echo "Inferred build file from 01g present"
	else
		echo "Problem: inferred build file from 01g is absent"
		exit 1
	fi

    # Main result prefixes
    local prefixes=(
        "hc_positive_control_untransformed"
        "no_correction_positive_control_untransformed"
        "pc_positive_control_untransformed"
    )

    local suffixes=(
        ".fastGWA.gz"
        "_manhattan.pdf"
        "_nocisChr_manhattan.pdf"
        "_qqplot.jpeg"
        "_nocisChr_qqplot.jpeg"
    )

    # Scatter plot method pairs
    local scatter_pairs=(
        "hc_no"
        "pc_hc"
        "pc_no"
    )

    while IFS=$'\t' read -r positive_control_cpg _; do
        [ -z "${positive_control_cpg}" ] && continue

        echo "Checking positive control: ${positive_control_cpg}"

        ## ---- Main outputs (per prefix) ----
        for prefix in "${prefixes[@]}"; do
            base="${section_01_dir}/01g/${prefix}_${positive_control_cpg}"

            for suffix in "${suffixes[@]}"; do
                file="${base}${suffix}"

                if [ ! -f "$file" ]; then
                    echo "Problem: missing file $file"
                    exit 1
                fi

                case "$suffix" in
                    ".fastGWA.gz")
                        echo "${prefix} ${positive_control_cpg} results present"
                        ;;
                    "_manhattan.pdf")
                        echo "${prefix} ${positive_control_cpg} Manhattan plot present"
                        ;;
                    "_nocisChr_manhattan.pdf")
                        echo "${prefix} ${positive_control_cpg} no cis chromosome Manhattan plot present"
                        ;;
                    "_qqplot.jpeg")
                        echo "${prefix} ${positive_control_cpg} QQ plot present"
                        ;;
                    "_nocisChr_qqplot.jpeg")
                        echo "${prefix} ${positive_control_cpg} no cis chromosome QQ plot present"
                        ;;
                esac
            done
        done

        ## ---- Scatter plots ----
        for pair in "${scatter_pairs[@]}"; do
            scatter_file="${section_01_dir}/01g/${positive_control_cpg}_${pair}_scatter.pdf"

            if [ ! -f "$scatter_file" ]; then
                echo "Problem: missing scatter plot $scatter_file"
                exit 1
            fi

            echo "scatter plot ${positive_control_cpg} ${pair} present"
        done

    done < <(
        awk -F'\t' '
            NR > 1 {
                gsub(/^[ \t]+|[ \t]+$/, "", $1)
                if ($1 != "" && $1 !~ /^#/)
                    print $1
            }
        ' "${methylation_processed_dir}/mqtl_pos_ctr_1g_filt.tsv"
    )
}

check_results_01 () {
	check_results_01a
	check_results_01b
	check_results_01c
	check_results_01d
	check_results_01e

    if [[ ! -z "${idat_directory}" ]]; then
        check_results_01f
    else
        echo "Skip check_results_01f: idat_directory not provided in config file"
    fi

    if [[ ! -z "${vcf_dir}" ]]; then
        check_results_01g
    else
        echo "Skip check_results_01g: vcf_dir not provided in config file"
    fi
    
}

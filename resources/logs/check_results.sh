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

	if [ -f "${snpchrtxt}" ] && [ -f "${snpchrpdf}" ]; then
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
		echo "Sex prediction plot file is absent."
		exit 1
	fi
}

check_results_01b () {

	if [ -f "${section_01_dir}/inferred_build.txt" ]; then
		echo "Inferred build file present"
	else
		echo "Problem: inferred build file is absent"
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
    echo "Problem: Sample missingness file is still in section_01_dir (${section_01_dir}), please delete it manually"
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

    if [ -f "${section_01_dir}/cohort_descriptives_commonids.RData" ]; then
		echo "cohort_descriptives_commonids.RData is present"
	else
		echo "Problem: cohort_descriptives_commonids.RData is absent"
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

	if [ -f "${smoking_pred_plot}" ]; then
		echo "Smoking prediction plot is present"
	else
		echo "Problem: Smoking prediction plot file not present"
		exit 1
	fi

	if [ -f "${section_01_dir}/age_prediction_correlation.png" ]; then
		echo "The matrix correlation plot among predicted age, age acceleration residuals, and chronological age is present"
	else
		echo "Problem: The matrix correlation plot of predicted age is absent"
		exit 1
	fi

	if ls "${section_01_dir}/pc_var_association_plot_${study_name}"*.pdf >/dev/null 2>&1; then
    	echo "The pc_var_association plot is present"
	else
    	echo "Problem: The pc_var_association plot is absent"
    	exit 1
	fi

	if ls "${section_01_dir}/pc_var_association_plot_${study_name}"*.pdf >/dev/null 2>&1; then
    	echo "The pc_var_association plot is present"
	else
    	echo "Problem: The pc_var_association plot is absent"
    	exit 1
	fi

	if [ -f "${section_01_dir}/edited_phenotypes_summary_${study_name}.Rdata" ] && [ -f "${section_01_dir}/edited_phenotype_distribution_plot_${study_name}.pdf" ]; then
		echo "Edited phenotypes summary Rdata file is present"
	else
		echo "Problem: Edited phenotypes summary Rdata file is absent"
		exit 1
	fi
}

check_results_01d () {

    if [ -f "${methylation_processed_dir}/mqtl_pos_ctr_filt.tsv" ]; then
        echo "filtered positive control file is present"
	else
		echo "Problem: filtered positive control file is absent"
		exit 1
    fi
	
	echo "Filtered positive control file is present: ${methylation_processed_dir}/mqtl_pos_ctr_filt.tsv"
    while IFS=$'\t' read -r positive_control_cpg _; do
        [ -z "${positive_control_cpg}" ] && continue

        echo "Checking positive control: ${positive_control_cpg}"

        if [ -f "${section_01_dir}/positive_control_untransformed_${positive_control_cpg}.fastGWAmlm.gz" ]; then
            echo "positive control ${positive_control_cpg} results present"
        else
            echo "Problem: positive control ${positive_control_cpg} results file not present: ${section_01_dir}/positive_control_untransformed_${positive_control_cpg}.fastGWAmlm.gz"
			exit 1
        fi

        if [ -f "${section_01_dir}/positive_control_untransformed_${positive_control_cpg}_manhattan.pdf" ]; then
            echo "positive control ${positive_control_cpg} Manhattan plot present"
        else
            echo "Problem: positive control ${positive_control_cpg} Manhattan plot file not present"
			exit 1
        fi

        if [ -f "${section_01_dir}/positive_control_untransformed_${positive_control_cpg}_nocisChr_manhattan.pdf" ]; then
            echo "positive control ${positive_control_cpg} no cis chromosome Manhattan plot present"
        else
            echo "Problem: positive control ${positive_control_cpg} no cis chromosome Manhattan plot file not present"
			exit 1
        fi

        if [ -f "${section_01_dir}/positive_control_untransformed_${positive_control_cpg}_qqplot.jpeg" ]; then
            echo "positive control ${positive_control_cpg} QQ plot present"
        else
            echo "Problem: positive control ${positive_control_cpg} QQ plot file not present"
			exit 1
        fi

        if [ -f "${section_01_dir}/positive_control_untransformed_${positive_control_cpg}_nocisChr_qqplot.jpeg" ]; then
            echo "positive control ${positive_control_cpg} no cis chromosome QQ plot present"
        else
            echo "Problem: positive control ${positive_control_cpg} no cis chromosome QQ plot file not present"
			exit 1
        fi

        # negative control corresponding name: NEG_${positive_control_cpg}
        neg="NEG_${positive_control_cpg}"
        echo "Checking negative control: ${neg}"

        if [ -f "${section_01_dir}/negative_control_untransformed_${neg}.fastGWAmlm.gz" ]; then
            echo "negative control ${neg} results present"
        else
            echo "Problem: negative control ${neg} results file not present: ${section_01_dir}/negative_control_untransformed_${neg}.fastGWAmlm.gz"
            exit 1
        fi

        if [ -f "${section_01_dir}/negative_control_untransformed_${neg}_manhattan.pdf" ]; then
            echo "negative control ${neg} Manhattan plot present"
        else
            echo "Problem: negative control ${neg} Manhattan plot file not present"
            exit 1
        fi

        if [ -f "${section_01_dir}/negative_control_untransformed_${neg}_nocisChr_manhattan.pdf" ]; then
            echo "negative control ${neg} no cis chromosome Manhattan plot present"
        else
            echo "Problem: negative control ${neg} no cis chromosome Manhattan plot file not present"
            exit 1
        fi

        if [ -f "${section_01_dir}/negative_control_untransformed_${neg}_qqplot.jpeg" ]; then
            echo "negative control ${neg} QQ plot present"
        else
            echo "Problem: negative control ${neg} QQ plot file not present"
            exit 1
        fi

        if [ -f "${section_01_dir}/negative_control_untransformed_${neg}_nocisChr_qqplot.jpeg" ]; then
            echo "negative control ${neg} no cis chromosome QQ plot present"
        else
            echo "Problem: negative control ${neg} no cis chromosome QQ plot file not present"
			exit 1
        fi

    done < <(awk -F'\t' 'NR>1 { gsub(/^[ \t]+|[ \t]+$/,"",$1); if($1!="" && $1!~ /^#/) print $1 }' "${methylation_processed_dir}/mqtl_pos_ctr_filt.tsv")

}

check_results_01e () {
    for i in {1..5}; do
        pc="PC${i}"

        if [ -f "${section_01_dir}/gwas_${pc}.fastGWAmlm.gz" ]; then
            echo "GWAS ${pc} results present"
        else
            echo "Problem: GWAS ${pc} results file not present: ${section_01_dir}/gwas_${pc}.fastGWAmlm.gz"
            exit 1
        fi

        if [ -f "${section_01_dir}/gwas_${pc}_manhattan_beta.pdf" ]; then
            echo "GWAS ${pc} beta Manhattan plot present"
        else
            echo "Problem: GWAS ${pc} beta Manhattan plot file not present: ${section_01_dir}/gwas_${pc}_manhattan_beta.pdf"
            exit 1
        fi
    done
}

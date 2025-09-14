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

	if [ -f "${section_01_dir}/easyQC_topmed.multi.AFCHECK.png" ]; then
		echo "easyQC plot present"
	else
		echo "Problem: easyQC plot is absent"
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

}

check_results_01c () {

	if [ -f "${section_01_dir}/methylation_summary.RData" ]; then
		echo "Methylation_summary.RData is present"
	else
		echo "Problem: methylation_summary.RData is absent"
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

	if [ "${measured_cellcounts}" != "NULL" ];then
		if [  -f "${section_01_dir}/cor_plot.pdf" ]; then
			echo "Correlation plot of observed vs predicted cell counts is present"
		else
			echo "Problem: correlation plot of observed vs predicted cell counts is absent"
			exit 1
        fi        

		if [  -f "${section_01_dir}/cor_matrix.txt" ]; then
			echo "Correlation matrix of observed vs predicted cell counts is present"
		else
			echo "Problem: correlation matrix of observed vs predicted cell counts is absent"
			exit 1
		fi
	else
		echo "Message: since measured_cellcounts are not provided, there is no output for cor_plot.pdf and cor_matrix.txt for observed vs predicted cell counts."

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

}

# ...existing code...

check_results_01d () {

    if [ ! -f "${methylation_processed_dir}/mqtl_pos_ctr_filt.tsv" ]; then
        echo "Problem: filtered positive control file is absent: ${methylation_processed_dir}/mqtl_pos_ctr_filt.tsv"
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

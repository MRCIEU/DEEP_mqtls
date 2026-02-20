#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

exec &> >(tee ${section_01e_logfile})
print_version

echo "Running GWAS for genetic PCs (PC1-PC10)"

covar_option=""
if [ -f "${ccovar_file}" ]; then
  ncols=$(awk -F'\t' 'NR==1{print NF; exit}' "${ccovar_file}" 2>/dev/null || echo 0)
  if [ "${ncols}" -gt 2 ]; then
    covar_option="--covar ${ccovar_file}"
  else
    echo "No categorical covariates detected in ${ccovar_file}; skipping --covar for GCTA."
  fi
else
  echo "Warning: ccovar file not found: ${ccovar_file}"
fi

if [ "${related}" = "yes" ]; then
    echo "Generating a new grm for mlm, to avoid proximal contamination"
    echo "Sorting variants using pgen"
    ${plink2} \
      --bfile ${bfile} \
      --make-pgen \
      --sort-vars \
      --new-id-max-allele-len 70 \
      --out ${bfile_sort} \
      --threads ${nthreads}

    ${plink2} \
      --pfile ${bfile_sort} \
      --new-id-max-allele-len 70 \
      --make-bed \
      --out ${bfile_sort}

    echo "Identifying SNPs that are highly LD with original PCA SNPs"
    ${plink} \
      --bfile ${bfile_sort} \
      --r2 gz \
      --ld-snp-list ${pca}.prune.in \
      --new-id-max-allele-len 70 \
      --ld-window-kb 1000 \
      --ld-window-r2 0.1 \
      --out "${high_ld}"

    echo "Creating a list of SNPs to exclude"
    #awk 'NR==1 && $6 ~ /[A-Za-z]/ {next} {print $6}' "${high_ld}.vcor" | awk 'NF' | sort -u > "${high_ld}.exclude.ids"
    zcat ${high_ld}.ld.gz |awk 'NR==1 && $6 ~ /[A-Za-z]/ {next} {print $6}' | awk 'NF' | sort -u > "${high_ld}.exclude.ids"

    echo "Excluding high LD regions in genome first"
    ${plink2} \
      --bfile ${bfile_sort} \
      --exclude range ${exclude_highld_region} \
      --new-id-max-allele-len 70 \
      --make-bed \
      --out ${bfile_sort}1

    mv ${bfile_sort}1.bed ${bfile_sort}.bed
		mv ${bfile_sort}1.bim ${bfile_sort}.bim
		mv ${bfile_sort}1.fam ${bfile_sort}.fam

    echo "Excluding highly-LD SNPs and pruning SNPs for GRM"
    ${plink2} \
      --bfile ${bfile_sort} \
	    --new-id-max-allele-len 70 \
	    --exclude "${high_ld}.exclude.ids" \
      --indep-pairwise 10000 5 0.1 \
	    --maf 0.2 \
	    --out ${snp_01e} \
	    --autosome \
	    --threads ${nthreads}

    echo "Generating new GRM for related individuals"
    ${plink2} \
      --bfile ${bfile} \
	    --new-id-max-allele-len 70 \
	    --extract "${snp_01e}.prune.in" \
	    --make-grm-bin \
	    --out "${grmfile_all_01e}" \
	    --threads "${nthreads}" \
	    --autosome

    echo "Sample size in creating kinship matrix: $(wc -l < "${grmfile_all_01e}.grm.id")"

    ${R_directory}Rscript resources/relateds/grm_distri.R \
	    "${grmfile_all_01e}" \
	    "${rel_cutoff}" \
	    "${grm_distribution}_01e"

    echo "Generating sparse GRM for fastGWA-mlm"
    ${gcta} \
	    --grm ${grmfile_all_01e} \
	    --make-bK-sparse 0.05 \
	    --autosome \
	    --make-grm \
	    --out ${grmfile_fast}_rel_01e \
	    --thread-num ${nthreads}
fi

for pc in {1..5}; do
    pc_col="PC${pc}"
    pheno_file="${home_directory}/processed_data/covariate_data/genetic_pc_gwas.PC${pc}.pheno"

    awk -v col="$((pc+2))" 'NR==1{next} {print $1, $2, $col}' ${genetic_pc_gwas} > ${pheno_file}

    echo "Running GWAS for ${pc_col}"
    #if [ "${related}" = "yes" ]; then
       # ${gcta} \
        #    --bfile "${bfile}" \
        #   --grm-sparse "${grmfile_fast}_rel_01e" \
        #    --fastGWA-mlm \
        #   --h2-limit 100 \
        #    --pheno "${pheno_file}" \
        #   --qcovar "${qcovar_noPC_file}" \
        #   ${covar_option} \
        #    --covar-maxlevel 300 \
        #    --out "${section_01_dir}/01e/gwas_${pc_col}" \
        #    --thread-num "${nthreads}"

    if [ "${related}" = "yes" ]; then
        ${gcta} \
            --bfile "${bfile}" \
            --mlma \
            --pheno "${pheno_file}" \
            --qcovar "${qcovar_noPC_file}" \
            ${covar_option} \
            --out "${section_01_dir}/01e/gwas_${pc_col}" \
            --thread-num "${nthreads}"
    
    tr -s " " < "${section_01_dir}/01e/gwas_${pc_col}.mlma" | gzip -c > "${section_01_dir}/01e/gwas_${pc_col}.mlma.gz"
    rm "${section_01_dir}/01e/gwas_${pc_col}.mlma"

    echo "make manhattan and qq plots for genetic ${pc_col}"
    echo "${section_01_dir}/01e/gwas_${pc_col}.mlma.gz" > "${section_01_dir}/01e/gwas_${pc_col}.file.txt"

    ${R_directory}Rscript resources/genetics/plot_gwas.R \
        "${section_01_dir}/01e/gwas_${pc_col}.file.txt" \
            9 \
            7 \
            1 \
            3 \
            2 \
            TRUE \
            0 \
            0 \
            0 \
            0 \
            beta


    elif [ "${related}" = "no" ]; then
        ${gcta} \
            --bfile "${bfile}" \
            --fastGWA-lr \
            --h2-limit 100 \
            --pheno "${pheno_file}" \
            --qcovar "${qcovar_noPC_file}" \
            ${covar_option} \
            --covar-maxlevel 300 \
            --out "${section_01_dir}/01e/gwas_${pc_col}" \
            --thread-num "${nthreads}"

    tr -s " " < "${section_01_dir}/01e/gwas_${pc_col}.fastGWA" | gzip -c > "${section_01_dir}/01e/gwas_${pc_col}.fastGWA.gz"
    rm "${section_01_dir}/01e/gwas_${pc_col}.fastGWA"

    echo "make manhattan and qq plots for genetic ${pc_col}"
    echo "${section_01_dir}/01e/gwas_${pc_col}.fastGWA.gz" > "${section_01_dir}/01e/gwas_${pc_col}.file.txt"

    ${R_directory}Rscript resources/genetics/plot_gwas.R \
        "${section_01_dir}/01e/gwas_${pc_col}.file.txt" \
            10 \
            8 \
            1 \
            3 \
            2 \
            TRUE \
            0 \
            0 \
            0 \
            0 \
            beta
    fi
done

echo "Successfully completed script 01e"
#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

exec &> >(tee ${section_01e_logfile})
print_version

echo "Running GWAS for genetic PCs (PC1-PC10)"

for pc in {1..2}; do
    pc_col="PC${pc}"
    pheno_file="genetic_pc_gwas.PC${pc}.pheno"

    awk -v col="$((pc+2))" 'NR==1{next} {print $1, $2, $col}' ${genetic_pc_gwas} > ${pheno_file}

    echo "Running GWAS for ${pc_col}"
    if [ "${related}" = "yes" ]; then
        ${gcta} \
            --bfile "${bfile}" \
            --grm-sparse "${grmfile_fast}_rel" \
            --fastGWA-mlm \
            --pheno "${pheno_file}" \
            --qcovar "${qcovar_noPC_file}" \
            --covar "${ccovar_file}" \
            --out "${section_01_dir}/gwas_${pc_col}" \
            --thread-num "${nthreads}"

    elif [ "${related}" = "no" ]; then
        ${gcta} \
            --bfile "${bfile}" \
            --grm-sparse "${grmfile_fast}_unrel" \
            --fastGWA-mlm \
            --pheno "${pheno_file}" \
            --qcovar "${qcovar_noPC_file}" \
            --covar "${ccovar_file}" \
            --out "${section_01_dir}/gwas_${pc_col}" \
            --thread-num "${nthreads}"
    fi

    tr -s " " < "${section_01_dir}/gwas_${pc_col}.fastGWA" | gzip -c > "${section_01_dir}/gwas_${pc_col}.fastGWAmlm.gz"
    rm "${section_01_dir}/gwas_${pc_col}.fastGWA"

    echo "make manhattan and qq plots for genetic PC ${pc_col}"
    echo "${section_01_dir}/gwas_${pc_col}.fastGWAmlm.gz" > "${section_01_dir}/gwas_${pc_col}.file.txt"

    ${R_directory}Rscript resources/genetics/plot_gwas.R \
        "${section_01_dir}/gwas_${pc_col}.file.txt" \
            0 \
            8 \
            1 \
            3 \
            2 \
            TRUE \
            0 \
            0 \
            0 \
            0

done

echo "Successfully completed script 01e"
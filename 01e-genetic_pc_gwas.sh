#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

exec &> >(tee ${section_01e_logfile})
print_version

echo "Please note this module is no longer required, since SNP loadings with LD-structure detection replace genetic PC GWAS as the default PCA LD-structure diagnostic."

echo "Running GWAS for genetic PCs (PC1-PC5)"

# Step 1: decide whether to use a sparse GRM and fastGWA mode; GRMs are generated in 01b
use_sparse_grm="no"
fastgwa_mode="lr"
grm_sparse_prefix=""

if [ "${related}" = "yes" ] && [ "${structured}" = "yes" ]; then
    echo "Related + Structured: using PC-Relate sparse GRM from 01b"
    fastgwa_mode="mlm"
    use_sparse_grm="yes"

    if [ -f "${grmfile_pcrelate}.grm.sp" ] && [ -f "${grmfile_pcrelate}.grm.id" ]; then
        grm_sparse_prefix="${grmfile_pcrelate}"
        echo "Using PC-Relate GCTA sparse GRM prefix: ${grm_sparse_prefix}"
    else
        echo "Error: Required PC-Relate sparse GRM files not found:"
        echo "  ${grmfile_pcrelate}.grm.sp"
        echo "  ${grmfile_pcrelate}.grm.id"
        exit 1
    fi

elif [ "${related}" = "yes" ] && [ "${structured}" = "no" ]; then
    echo "Related + non-structured: using KING sparse GRM from 01b"
    fastgwa_mode="mlm"
    use_sparse_grm="yes"
    grm_sparse_prefix="${grmfile_king}_sparse"

    if [ -f "${grmfile_king}_sparse.grm.sp" ] && [ -f "${grmfile_king}_sparse.grm.id" ]; then
        echo "Using KING sparse GRM prefix: ${grm_sparse_prefix}"
    else
        echo "Error: Required KING sparse GRM files not found:"
        echo "  ${grmfile_king}_sparse.grm.sp"
        echo "  ${grmfile_king}_sparse.grm.id"
        exit 1
    fi

elif [ "${related}" = "no" ]; then
    echo "Unrelated sample: no GRM needed (fastGWA-lr)"
    fastgwa_mode="lr"
    use_sparse_grm="no"
    grm_sparse_prefix=""

else
    echo "Error: Invalid related/structured flag combination. Set both to yes or no."
    exit 1
fi

fastgwa_option=()
grm_sparse_option=()
if [ "${use_sparse_grm}" = "yes" ]; then
    fastgwa_option=("--fastGWA-mlm")
    grm_sparse_option=("--grm-sparse" "${grm_sparse_prefix}")
elif [ "${fastgwa_mode}" = "lr" ]; then
    fastgwa_option=("--fastGWA-lr")
else
    echo "Error: Unsupported fastGWA mode: ${fastgwa_mode}"
    exit 1
fi

covar_args=()
if [ -f "${ccovar_file}" ]; then
  ncols=$(awk -F'\t' 'NR==1{print NF; exit}' "${ccovar_file}" 2>/dev/null || echo 0)
  if [ "${ncols}" -gt 2 ]; then
    covar_args=("--covar" "${ccovar_file}")
  else
    echo "No categorical covariates detected in ${ccovar_file}; skipping --covar for GCTA."
  fi
else
  echo "Warning: ccovar file not found: ${ccovar_file}"
fi

for pc in {1..5}; do
    pc_col="PC${pc}"
    pheno_file="${home_directory}/processed_data/covariate_data/genetic_pc_gwas.PC${pc}.pheno"

    awk -v col="$((pc+2))" 'NR==1{next} {print $1, $2, $col}' ${genetic_pc_gwas} > ${pheno_file}

    echo "Running GWAS for ${pc_col}"
    if [ "${use_sparse_grm}" = "yes" ]; then
      echo "Using sparse GRM prefix: ${grm_sparse_prefix}"
    else
      echo "No GRM mode: fastGWA-lr"
    fi

    ${gcta} \
      --bfile "${bfile}" \
      "${grm_sparse_option[@]}" \
      "${fastgwa_option[@]}" \
      --h2-limit 100 \
      --pheno "${pheno_file}" \
      --qcovar "${qcovar_noPC_file}" \
      "${covar_args[@]}" \
      --covar-maxlevel 500 \
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

done

echo "Successfully completed script 01e"
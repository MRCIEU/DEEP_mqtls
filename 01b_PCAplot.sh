#!/bin/bash
# DEEP global PCA v2 for the DEEP profile paper (branch global_pca_v2).
#
# Projects the cohort onto a reference PCA (HGDP + 1000 Genomes) built only on SNPs that
# all participating DEEP cohorts have, so every cohort has full coverage.
#
#   ./01b_PCAplot.sh -c <config>
#   ./check_upload.sh -c <config> 01b_pca upload        (or: 01b_pca upload manual)
#
# Optional config entries (see resources/global_pca/parameters_01b_pca):
#   pca_mode="post01b"          default; for cohorts that have completed 01b
#   pca_mode="genetic_only"     for cohorts without methylation data yet; needs bfile_raw,
#                               quality_scores and pca_sample_list (FID IID of the samples
#                               selected for methylation profiling)
#   pca_share_level="individual"  default; uploads PC1-PC10 per sample with no IDs, rows shuffled
#   pca_share_level="aggregate"   uploads only sample counts in fixed bins (counts < 5 hidden)
#
# Results go to results/01b_pca (uploaded). Scores with IDs stay in
# processed_data/genetic_data/global_pca_v2 and are never uploaded.

source resources/setup.sh "$@"
set -- $concatenated
source resources/global_pca/parameters_01b_pca

mkdir -p "${section_01b_pca_dir}/logs" "${gpca_work_dir}"
exec &> >(tee "${section_01b_pca_logfile}")
print_version
echo "01b_PCAplot: DEEP global PCA v2"
echo "pca_mode=${pca_mode}  pca_share_level=${pca_share_level}"

if [ "${pca_mode}" != "post01b" ] && [ "${pca_mode}" != "genetic_only" ]; then
    echo "Error: pca_mode must be 'post01b' or 'genetic_only' (config: ${config_file})" >&2
    exit 1
fi
if [ "${pca_share_level}" != "individual" ] && [ "${pca_share_level}" != "aggregate" ]; then
    echo "Error: pca_share_level must be 'individual' or 'aggregate' (config: ${config_file})" >&2
    exit 1
fi
for f in "${gpca_ref_prefix}_loadings.tsv.bgz" "${gpca_ref_prefix}_scores_flat.tsv.gz" \
         "${gpca_ref_prefix}_bins.tsv" "${gpca_ref_prefix}.snplist"; do
    if [ ! -s "${gpca_ref_dir}/${f}" ]; then
        echo "Error: reference file missing: ${gpca_ref_dir}/${f}. Run 'git pull' on the global_pca_v2 branch." >&2
        exit 1
    fi
done

# reference SNP positions (1-based ranges) to subset genotypes early
awk -F'[:_]' -v OFS='\t' '{print $1, $2, $2, $0}' "${gpca_ref_dir}/${gpca_ref_prefix}.snplist" > "${gpca_ref_ranges}"
echo "Reference SNPs: $(wc -l < "${gpca_ref_ranges}")"

if [ "${pca_mode}" = "post01b" ]; then
    # =================================================================
    # post01b: cleaned 01b genotypes (GRCh38, chr:pos_A1_A2 IDs)
    # =================================================================
    for ext in bed bim fam; do
        if [ ! -s "${bfile}.${ext}" ]; then
            echo "Error: ${bfile}.${ext} not found. Run 01b first, or set pca_mode=\"genetic_only\" in the config." >&2
            exit 1
        fi
    done
    id_ok=$(head -n 2000 "${bfile}.bim" | awk '$2 ~ /^[0-9]+:[0-9]+_[ACGT]+_[ACGT]+$/ {n++} END {print (NR > 0 && n / NR > 0.9) ? "yes" : "no"}')
    if [ "${id_ok}" != "yes" ]; then
        echo "Error: SNP IDs in ${bfile}.bim are not chr:pos_A1_A2; this does not look like 01b output." >&2
        exit 1
    fi
    build="38"
    ${plink2} \
        --bfile "${bfile}" \
        --extract bed1 "${gpca_ref_ranges}" \
        --new-id-max-allele-len 500 \
        --make-bed \
        --output-chr 26 \
        --out "${gpca_bfile}" \
        --threads ${nthreads}

else
    # =================================================================
    # genetic_only: same SNP QC as 01b, for the samples in pca_sample_list
    # =================================================================
    if [ ! -s "${pca_sample_list}" ]; then
        echo "Error: pca_sample_list not found or empty: ${pca_sample_list}" >&2
        echo "Provide FID IID (one sample per line) of the samples selected for methylation profiling." >&2
        exit 1
    fi
    for ext in bed bim fam; do
        if [ ! -s "${bfile_raw}.${ext}" ]; then
            echo "Error: ${bfile_raw}.${ext} not found (bfile_raw in config)." >&2
            exit 1
        fi
    done
    if [ ! -s "${quality_scores}" ]; then
        echo "Error: quality_scores file not found: ${quality_scores}" >&2
        exit 1
    fi
    n_list=$(wc -l < "${pca_sample_list}")
    # two columns: match FID and IID; one column: match IID (as plink2 --keep does)
    n_match=$(awk 'NR == FNR {if (NF >= 2) k[$1" "$2] = 1; else if (NF == 1) k[$1] = 1; next}
                   (($1" "$2) in k) || ($2 in k) {n++} END {print n + 0}' \
        "${pca_sample_list}" "${bfile_raw}.fam")
    echo "pca_sample_list: ${n_list} samples, ${n_match} found in ${bfile_raw}.fam"
    if [ "${n_match}" -eq 0 ]; then
        echo "Error: none of the samples in ${pca_sample_list} are in ${bfile_raw}.fam." >&2
        exit 1
    fi

    echo "Removing SNPs with INFO < 0.80 before liftover"
    awk '$3 < 0.80 {print $1}' "${quality_scores}" > "${gpca_bfile}.lowinfoSNPs.txt"
    ${plink2} --bfile "${bfile_raw}" \
        --keep "${pca_sample_list}" \
        --exclude "${gpca_bfile}.lowinfoSNPs.txt" \
        --new-id-max-allele-len 500 \
        --make-bed \
        --output-chr 26 \
        --out "${gpca_bfile}_info_filtered" \
        --threads ${nthreads}

    echo "Inferring genome build and running liftover if necessary (config genome_build: ${genome_build})"
    rm -f "${gpca_miss_liftover}"
    ${R_directory}Rscript resources/datacheck/liftover.R \
        "${gpca_bfile}_info_filtered" \
        ${genome_build} \
        "${gpca_miss_liftover}" \
        "${gpca_liftover_map}" \
        "${gpca_work_dir}" \
        "01b_pca"
    build=$(cat "${gpca_work_dir}/01b_pca_inferred_build.txt")
    echo "Inferred build: ${build}"

    if [ "${build}" -eq 37 ]; then
        exclude_args=()
        if [ -f "${gpca_miss_liftover}" ]; then
            echo "Excluding SNPs that cannot be lifted over"
            exclude_args=(--exclude "${gpca_miss_liftover}")
        fi
        ${plink2} --bfile "${gpca_bfile}_info_filtered" \
            "${exclude_args[@]}" \
            --update-map "${gpca_liftover_map}" \
            --new-id-max-allele-len 500 \
            --make-bed \
            --output-chr 26 \
            --out "${gpca_bfile}_b38" \
            --threads ${nthreads}
        ${plink2} --bfile "${gpca_bfile}_b38" --sort-vars --make-pgen --out "${gpca_bfile}_b38_sorted" --threads ${nthreads}
        ${plink2} --pfile "${gpca_bfile}_b38_sorted" --make-bed --output-chr 26 --out "${gpca_bfile}_b38" --threads ${nthreads}
        rm -f "${gpca_bfile}_b38_sorted".{pgen,psam,pvar,log}
    elif [ "${build}" -eq 38 ]; then
        ${plink2} --bfile "${gpca_bfile}_info_filtered" \
            --new-id-max-allele-len 500 \
            --make-bed \
            --output-chr 26 \
            --out "${gpca_bfile}_b38" \
            --threads ${nthreads}
    else
        echo "Error: could not infer genome build (got '${build}')." >&2
        exit 1
    fi
    rm -f "${gpca_bfile}_info_filtered".{bed,bim,fam}

    echo "SNP QC on reference SNP positions (MAF ${snp_maf}, HWE ${snp_hwe}, SNP missingness ${snp_miss}, sample missingness ${snp_imiss})"
    ${plink2} \
        --bfile "${gpca_bfile}_b38" \
        --extract bed1 "${gpca_ref_ranges}" \
        --maf ${snp_maf} \
        --hwe ${snp_hwe} \
        --geno ${snp_miss} \
        --mind ${snp_imiss} \
        --chr 1-22 \
        --new-id-max-allele-len 500 \
        --make-bed \
        --output-chr 26 \
        --out "${gpca_bfile}_format" \
        --threads ${nthreads}
    rm -f "${gpca_bfile}_b38".{bed,bim,fam}

    echo "Formatting SNP IDs to chr:pos_A1_A2"
    ${plink2} \
        --bfile "${gpca_bfile}_format" \
        --set-all-var-ids @:#_\$1_\$2 \
        --new-id-max-allele-len 500 \
        --make-bed \
        --output-chr 26 \
        --out "${gpca_bfile}_ids" \
        --threads ${nthreads}

    rm -f "${gpca_snpfail}"; touch "${gpca_snpfail}"
    ${R_directory}Rscript resources/genetics/harmonization.R "${gpca_bfile}_ids.bim" "${gpca_snpfail}"
    sort -u "${gpca_snpfail}" > "${gpca_snpfail}.sorted"

    ${plink2} \
        --bfile "${gpca_bfile}_ids" \
        --exclude "${gpca_snpfail}.sorted" \
        --rm-dup exclude-mismatch \
        --new-id-max-allele-len 500 \
        --make-bed \
        --output-chr 26 \
        --out "${gpca_bfile}" \
        --threads ${nthreads}
    rm -f "${gpca_bfile}_format".{bed,bim,fam} "${gpca_bfile}_ids".{bed,bim,fam}
fi

echo "Genotypes for projection: $(wc -l < "${gpca_bfile}.bim") SNPs, $(wc -l < "${gpca_bfile}.fam") samples"

version=$(git describe --tags --always 2>/dev/null || echo unknown)
echo "Projecting onto the DEEP global PCA v2 reference"
${Python_directory}python "${scripts_directory}/resources/global_pca/project_global_pca.py" \
    --bfile "${gpca_bfile}" \
    --ref-dir "${gpca_ref_dir}" \
    --ref-prefix "${gpca_ref_prefix}" \
    --study "${study_name}" \
    --mode "${pca_mode}" \
    --share-level "${pca_share_level}" \
    --build "${build}" \
    --out-dir "${section_01b_pca_dir}" \
    --work-dir "${gpca_work_dir}" \
    --pipeline-version "${version}" \
    --threads ${nthreads} \
    --mem-gb ${mem}

echo "Checking results/01b_pca for unexpected files and sample IDs"
id_files=("${gpca_work_dir}/projected_ids.txt" "${gpca_bfile}.fam")
if [ "${pca_mode}" = "genetic_only" ]; then
    id_files+=("${pca_sample_list}")
else
    id_files+=("${bfile}.fam")
fi
${Python_directory}python "${scripts_directory}/resources/global_pca/check_01b_pca_outputs.py" \
    --out-dir "${section_01b_pca_dir}" \
    --study "${study_name}" \
    --share-level "${pca_share_level}" \
    --id-files "${id_files[@]}"

echo ""
echo "Results: ${section_01b_pca_dir}"
echo "Please look at ${section_01b_pca_dir}/${study_name}_globalPCA_v2.png before uploading."
echo "Upload with: ./check_upload.sh -c ${config_file} 01b_pca upload"
echo "Successfully completed script 01b_PCAplot"

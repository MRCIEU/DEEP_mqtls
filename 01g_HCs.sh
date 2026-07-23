#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

exec &> >(tee ${section_01g_logfile})
print_version

containsElement () {
  local e
  for e in "${@:2}"; do [[ "$e" == "$1" ]] && return 0; done
  echo "There is no method for ${1}."
  echo "Please run:"
  echo "./01g_HCs.sh [arg]"
  echo "where arg is an optional argument that can be one of:"
  printf '%s\n' ${@:2}
  return 1
}

arg="all"
declare -a sections=('all' 'vcf' 'hc' 'gwas')

if [ -n "${1}" ]; then
  arg="${1}"
  containsElement "${1}" ${sections[@]}
fi

section_message () {
  echo "-----------------------------------------------"
  echo ""
  echo "$1 section"
  echo ""
  echo "to run this part on its own type:"
  echo "$ ./01g_HCs.sh $1"
  echo ""
  echo "-----------------------------------------------"
  echo ""
  echo ""
}

# because PCs are constructed from total variance explained, many are associated with genetically localized structures rather than population structure.
# computing haplotype components (HCs) through PBWTpaint, which captures fine-scale local haplotype sharing patterns, can better account for subtle population structure in GWAS, but not for relatedness.

if [ "$arg" = "vcf" ] || [ "$arg" = "all" ]
then
section_message "vcf"

echo "Checking VCF files in: ${vcf_dir}"
# exclude Chr X
missing=()
for i in {1..22}; do
  file=$(find "${vcf_dir}" -maxdepth 1 -type f \
    -regextype posix-extended \
    -iregex ".*/chr${i}([^0-9].*)?\.vcf(\.(gz|bgz|zip))?" \
    -print -quit)

  if [ -z "${file}" ]; then
    echo "MISSING: chr${i}"
    missing+=("chr${i}")
  else
    echo "FOUND:  chr${i} -> ${file}"
  fi
done

if [ ${#missing[@]} -gt 0 ]; then
  echo "ERROR: missing chromosomes: ${missing[*]}" >&2
  exit 2
else
  echo "All chr1..chr22 VCF files present."
fi

for i in {1..22}; do

    chr_vcf=$(find "${vcf_dir}" -maxdepth 1 -type f \
    -regextype posix-extended \
    -iregex ".*/chr${i}([^0-9].*)?\.vcf(\.(gz|bgz|zip))?" \
    -print -quit)

  if [ -z "${chr_vcf}" ]; then
    echo "ERROR: could not locate VCF for chr${i} in ${vcf_dir}"
    exit 2
  fi

  echo "Processing chr${i} VCF: ${chr_vcf}"
  echo "Converting to PLINK bed/bim/fam format for liftover check"

  ${plink2} --vcf "${chr_vcf}" \
    --new-id-max-allele-len 500 \
    --set-all-var-ids @:#_\$1_\$2 \
    --make-bed \
    --out "${genetic_processed_dir}/chr${i}_tmp"

  # per-chromosome liftover temp files
  temp_miss_liftover_chr="${section_01_dir}/01g_chr${i}_miss_liftover.txt"
  temp_liftover_map_chr="${section_01_dir}/01g_chr${i}_liftover_map.txt"

  echo "Confirming genome build and performing liftover"
  ${R_directory}Rscript resources/datacheck/liftover.R \
    "${genetic_processed_dir}/chr${i}_tmp" \
    "${genome_build}" \
    "${temp_miss_liftover_chr}" \
    "${temp_liftover_map_chr}" \
    "${section_01_dir}" \
    "01g"

  inferred_build=$(cat "${section_01_dir}/01g_inferred_build.txt")
  output_vcf="${genetic_processed_dir}/chr${i}_lifted"

  if [ "${inferred_build}" -eq 37 ]; then
    echo "Liftover chr${i} (build 37 -> 38)"
    if [ -f "${temp_miss_liftover_chr}" ]; then
      ${plink2} --vcf "${chr_vcf}" \
        --new-id-max-allele-len 500 \
        --set-all-var-ids @:#_\$1_\$2 \
        --exclude "${temp_miss_liftover_chr}" \
        --update-map "${temp_liftover_map_chr}" \
        --export vcf bgz \
        --out "${output_vcf}"
    else
      ${plink2} --vcf "${chr_vcf}" \
        --new-id-max-allele-len 500 \
        --set-all-var-ids @:#_\$1_\$2 \
        --update-map "${temp_liftover_map_chr}" \
        --export vcf bgz \
        --out "${output_vcf}"
    fi
  elif [ "${inferred_build}" -eq 38 ]; then
    echo "Build 38 detected; exporting vcf with rsid formatted"
    ${plink2} --vcf "${chr_vcf}" \
      --new-id-max-allele-len 500 \
      --set-all-var-ids @:#_\$1_\$2 \
      --export vcf bgz \
      --out "${output_vcf}"
  else
    echo "ERROR: Unsupported inferred build '${inferred_build}'"; exit 1
  fi

  # clean up temp files
  rm "${genetic_processed_dir}/chr${i}_tmp.bed"
  rm "${genetic_processed_dir}/chr${i}_tmp.bim"
  rm "${genetic_processed_dir}/chr${i}_tmp.fam"
  rm "${genetic_processed_dir}/chr${i}_tmp.log"

done

# then vcf files needs to be cleaned with SNPs with samples and variants from 01b
# input vcf from process_data folder
# input sample and snp list from processed data of 01b
# to confirm snp format here
echo "Preparing sample and SNP lists from 01b processed data"
awk '{print $2}' "${bfile}.bim" > "${genetic_processed_dir}/snp_list"
awk '{print $2}' "${bfile}.fam" > "${genetic_processed_dir}/sample_list"

echo "$(wc -l < "${genetic_processed_dir}/snp_list") SNPs in snp_list"
echo "$(wc -l < "${genetic_processed_dir}/sample_list") samples in sample_list"

for i in {1..22}; do
    echo "Preparing cleaned chr${i} vcf data"
        ${plink2} --vcf "${genetic_processed_dir}/chr${i}_lifted.vcf.gz" \
            --keep "${genetic_processed_dir}/sample_list" \
            --extract "${genetic_processed_dir}/snp_list" \
            --export vcf bgz \
            --out "${genetic_processed_dir}/chr${i}_data"
done

for i in {1..22}; do
    rm "${genetic_processed_dir}/chr${i}_lifted.vcf.gz"
done

echo "Successfully completed script 01g vcf chunk"
fi

if [ "$arg" = "hc" ] || [ "$arg" = "all" ]
then
section_message "hc"

for i in {1..22}; do

    echo "Computing pair-wise genetic distance for chr${i}..."

    ${pbwt} -readVcfGT "${genetic_processed_dir}/chr${i}_data.vcf.gz" \
            -paintSparse "${genetic_processed_dir}/chr${i}_data" 100 2 0

done

cp "${scripts_directory}/resources/genetics/combine_chunklength.cpp" "${home_directory}/processed_data/genetic_data/combine_chunklength_edit.cpp"
cp "${scripts_directory}/resources/genetics/gzstream.C" "${home_directory}/processed_data/genetic_data/gzstream.C"
cp "${scripts_directory}/resources/genetics/gzstream.h" "${home_directory}/processed_data/genetic_data/gzstream.h"

# combine_chunklength_edit.cpp
# needs number of SNPs on each chromosome

counts=()

for i in {1..22}; do
    vcf="${genetic_processed_dir}/chr${i}_data.vcf.gz"
    out="${genetic_processed_dir}/chr${i}_temp_plink"

    ${plink2} --vcf "$vcf" --freq --out "$out" --threads "${nthreads}" >/dev/null 2>&1

    n=$(($(wc -l < "${out}.afreq") - 1))

    counts+=("$n")
    rm -f "${out}.afreq"
    rm -f "${out}.log"
done

counts_list=$(printf "%s, " "${counts[@]}")
counts_list=${counts_list%, } 
echo "comma-separated list:"
echo "$counts_list"

# needs number of individuals
sample_size=$(wc -l < "${genetic_processed_dir}/sample_list")

echo "Sample size (number of individuals): ${sample_size}"

# needs path of first chunklength file

echo "Modifying combine_chunklength_edit.cpp with counts, sample size and its first file"

sed -i "s|vector<double> total_SNP = {};|vector<double> total_SNP = { ${counts_list} };|" \
  "${home_directory}/processed_data/genetic_data/combine_chunklength_edit.cpp"

sed -i "s/int nind=N;/int nind=${sample_size};/" \
"${home_directory}/processed_data/genetic_data/combine_chunklength_edit.cpp"

sed -i "s|PROCESS_DIR|${genetic_processed_dir}|g" "${home_directory}/processed_data/genetic_data/combine_chunklength_edit.cpp"

sed -i 's|ogzstream outFile("");|ogzstream outFile("'"${home_directory}/processed_data/genetic_data/full_chunklength.txt.gz"'");|' \
  "${home_directory}/processed_data/genetic_data/combine_chunklength_edit.cpp"

echo "$RUN_CMD g++ compiling combine_chunklength_edit.cpp"

$RUN_CMD bash -c "
g++ ${home_directory}/processed_data/genetic_data/combine_chunklength_edit.cpp \
    -o ${home_directory}/processed_data/genetic_data/combine \
    -I\$CONDA_PREFIX/include \
    -L\$CONDA_PREFIX/lib \
    -lz -lpthread -llapack -lblas -std=c++0x -g -O3
"

echo "Combining chunklength files to generate full chunklength file"
${home_directory}/processed_data/genetic_data/combine

# default generate top 100 HCs
# if sample size < 100, then number of HCs will be the sample size; variance explained plot
echo "Generating haplotype components (HCs)"
${R_directory}Rscript resources/genetics/generate_hcs.R \
    "${full_chunklength}" \
    "${bfile}.fam" \
    "${hcs_file}" \
    "${hcs_kernel_prop}" \
    "${hcs_plot}" \
    "${hcs_scatter_plot}"

echo "Selecting HCs using profile likelihood"
${R_directory}Rscript resources/genetics/select_hcs_profile_likelihood.R \
    "${hcs_kernel_prop}" \
    "${hcs_file}" \
    "${hcs_profile_selected_file}" \
    "${hcs_profile_likelihood}" \
    "${hcs_profile_selection}" \
    "${hcs_profile_plot}"

echo "Generating side-by-side scatter plots for HC1vHC2, HC2vHC3 and HC1vHC3, coloured by factors (population-structure check)"
${R_directory}Rscript resources/genetics/hc_pop_scatter.R \
    "${hcs_file}" \
    "${covariates_combined}.txt" \
    "${winsorized_phenotype_file}" \
    "${study_name}" \
    "${section_01_dir}"


echo "Computing correlation and linear model between HCs and PCs"
# correlation betweene HCs and PCs; linear reg for hcs against pcs; pcs against hcs
${R_directory}Rscript resources/genetics/corr_lm_hc_pc.R \
    ${pcs_all} \
    ${hcs_file} \
    ${genetic_outlier_ids} \
    ${section_01_dir}

# Generate qcovar file with profile-likelihood selected HCs for GWAS.
${R_directory}Rscript resources/methylation/generate_qcovar_with_hcs.R \
    "${covariates_combined}.txt" \
    "${bfile}.fam" \
    "${hcs_profile_selected_file}" \
    "${qcovar_hc_file}" \
    "${scripts_directory}" 

echo "Successfully completed script 01g hc chunk"
fi

# qcovar file with PCs are from 01d
# qcovar file with no correction from  01e

if [ "$arg" = "gwas" ] || [ "$arg" = "all" ]
then
section_message "gwas"

echo "Checking positive control file for 01g"
${R_directory}Rscript resources/genetics/check_positive_controls.R \
    "${positive_control_1g_file}" \
    "${filt_positive_control_1g_file}" \
    "${com_id}" \
    "${methylation_no_outliers_gwas}"

# CpG for pop structure
# for hcs corrected gwas results, generate figures # gwas manhattan and qq
# phenotype sensitive to population structure
# epistructure
# perform GWAS
echo "Performing GWAS with HCs/PCs/no-correction as covariates, other covs, like cell types, age and sex still included"
# GWAS with HCs as covariates

# choose sparse GRM and fastGWA mode the same way as 01d/01e
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

if [ "${#covar_args[@]}" -gt 0 ]; then
  echo "Using categorical covariates from ${ccovar_file}"
else
  echo "No categorical covariates passed to GCTA"
fi

base_methylation_no_outliers="${methylation_no_outliers_gwas%.Robj}"

tail -n +2 "${filt_positive_control_1g_file}" | while IFS=$'\t' read -r positive_control_cpg positive_control_snp_chr positive_control_snp_pos rsid positive_control_snp_window positive_control_threshold
do
    echo "Processing positive control: $positive_control_cpg, SNP chr: $positive_control_snp_chr, pos: $positive_control_snp_pos, rsid: $rsid, window: $positive_control_snp_window, threshold: $positive_control_threshold"

    echo "Extracting methylation values for positive control CpG ${positive_control_cpg}"

    ${R_directory}Rscript -e "load('${methylation_no_outliers_gwas}'); \
                if (!exists('norm.beta')) norm.beta <- get(ls()[1]); \
                row <- norm.beta[rownames(norm.beta) == '${positive_control_cpg}', , drop=FALSE]; \
                ids <- colnames(norm.beta); \
                vals <- as.numeric(row[1,]); \
                df <- data.frame(IID = ids, stringsAsFactors = FALSE); \
                df[['${positive_control_cpg}']] <- vals; \
                write.table(df, file='${base_methylation_no_outliers}.${positive_control_cpg}.positive_control', sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)"

    echo "making gcta input for positive control CpG ${positive_control_cpg} (untransformed)"
    nrow=`cat ${base_methylation_no_outliers}.${positive_control_cpg}.positive_control | wc -l`
    if [ "$nrow" -lt "2" ]; then
        echo "The positive control CpG ${positive_control_cpg} appears to be missing for mQTL analysis on untransformed methylation data. Please check."
    else 
        ${R_directory}Rscript resources/genetics/make_control.R \
            "${base_methylation_no_outliers}.${positive_control_cpg}.positive_control" \
            "${intersect_ids_plink}" \
            "${base_methylation_no_outliers}.${positive_control_cpg}.positive_control.gcta"

        echo "Perform fastGWA (untransformed) in positive control"
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
          --pheno "${base_methylation_no_outliers}.${positive_control_cpg}.positive_control.gcta" \
          --qcovar "${qcovar_file}" \
          "${covar_args[@]}" \
          --covar-maxlevel 300 \
          --out "${section_01_dir}/01g/pc_positive_control_untransformed_${positive_control_cpg}" \
          --thread-num "${nthreads}"

        ${gcta} \
          --bfile "${bfile}" \
          "${grm_sparse_option[@]}" \
          "${fastgwa_option[@]}" \
          --h2-limit 100 \
          --pheno "${base_methylation_no_outliers}.${positive_control_cpg}.positive_control.gcta" \
          --qcovar "${qcovar_hc_file}" \
          "${covar_args[@]}" \
          --covar-maxlevel 300 \
          --out "${section_01_dir}/01g/hc_positive_control_untransformed_${positive_control_cpg}" \
          --thread-num "${nthreads}"

        ${gcta} \
          --bfile "${bfile}" \
          "${grm_sparse_option[@]}" \
          "${fastgwa_option[@]}" \
          --h2-limit 100 \
          --pheno "${base_methylation_no_outliers}.${positive_control_cpg}.positive_control.gcta" \
          --qcovar "${qcovar_noPC_file}" \
          "${covar_args[@]}" \
          --covar-maxlevel 300 \
          --out "${section_01_dir}/01g/no_correction_positive_control_untransformed_${positive_control_cpg}" \
          --thread-num "${nthreads}"

        for prefix in pc_positive_control_untransformed_${positive_control_cpg} \
                      hc_positive_control_untransformed_${positive_control_cpg} \
                      no_correction_positive_control_untransformed_${positive_control_cpg}; do

            in_txt="${section_01_dir}/01g/${prefix}.fastGWA"
            out_gz="${section_01_dir}/01g/${prefix}.fastGWA.gz"

            if [ -f "${in_txt}" ]; then
                tr -s " " < "${in_txt}" | gzip -c > "${out_gz}"
                rm "${in_txt}"
            else
                echo "Warning: expected fastGWA output not found: ${in_txt}"
            fi
        done

        echo "make manhattan and qq plots"
        {
          echo "${section_01_dir}/01g/pc_positive_control_untransformed_${positive_control_cpg}.fastGWA.gz"
          echo "${section_01_dir}/01g/hc_positive_control_untransformed_${positive_control_cpg}.fastGWA.gz"
          echo "${section_01_dir}/01g/no_correction_positive_control_untransformed_${positive_control_cpg}.fastGWA.gz"
        } > "${section_01_dir}/01g/positive.control.untransformed.file.txt"

        ${R_directory}Rscript resources/genetics/plot_gwas.R \
            "${section_01_dir}/01g/positive.control.untransformed.file.txt" \
            10 \
            8 \
            1 \
            3 \
            2 \
            TRUE \
            "${positive_control_snp_chr}" \
            "${positive_control_snp_pos}" \
            "${positive_control_snp_window}" \
            "${positive_control_threshold}" \
            p

        echo "Comparing GWAS results across PCs, HCs and non-corrected results by scatter plots"
        ${R_directory}Rscript resources/genetics/comp_gwas_results.R \
          "${section_01_dir}/01g/positive.control.untransformed.file.txt" \
          "${positive_control_cpg}" \
          10 \
          8 \
          9 \
          1 \
          3 \
          2 \
          TRUE \
          "${positive_control_snp_chr}" \
          "${positive_control_snp_pos}" \
          "${positive_control_snp_window}" \
          "${positive_control_threshold}"

    fi

done

echo "Successfully completed script 01g gwas chunk"
fi

if [ "$arg" = "all" ]; then
  echo "Successfully completed script 01g"
fi

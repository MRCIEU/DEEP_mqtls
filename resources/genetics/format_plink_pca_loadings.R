args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(
    "Usage: format_plink_pca_loadings.R ",
    "<plink.eigenvec.var> <pca.prune.in> <output.tsv.gz> <n_pcs>"
  )
}

input_file <- args[1]
expected_snp_file <- args[2]
output_file <- args[3]
n_pcs <- as.integer(args[4])

if (!file.exists(input_file)) {
  stop("PLINK SNP-loading file does not exist: ", input_file)
}
if (!file.exists(expected_snp_file)) {
  stop("Expected PCA SNP file does not exist: ", expected_snp_file)
}
if (is.na(n_pcs) || n_pcs < 1L) {
  stop("n_pcs must be a positive integer.")
}

script_argument <- grep(
  "^--file=",
  commandArgs(trailingOnly = FALSE),
  value = TRUE
)
if (length(script_argument) != 1L) {
  stop("Could not determine the formatter script location.")
}
script_file <- sub("^--file=", "", script_argument)
source(file.path(dirname(normalizePath(script_file)), "pca_loadings.R"))

native_loadings <- utils::read.delim(
  input_file,
  header = TRUE,
  check.names = FALSE,
  comment.char = "",
  stringsAsFactors = FALSE
)

chromosome_column <- if ("#CHROM" %in% colnames(native_loadings)) {
  "#CHROM"
} else {
  "CHROM"
}
required_columns <- c(
  chromosome_column, "POS", "ID", "MAJ", "NONMAJ",
  paste0("PC", seq_len(n_pcs))
)
missing_columns <- setdiff(required_columns, colnames(native_loadings))
if (length(missing_columns) > 0L) {
  stop(
    "Missing columns in PLINK SNP-loading file: ",
    paste(missing_columns, collapse = ", ")
  )
}

output <- data.frame(
  CHR = native_loadings[[chromosome_column]],
  POS = native_loadings[["POS"]],
  SNP = as.character(native_loadings[["ID"]]),
  LOADING_ALLELE = as.character(native_loadings[["MAJ"]]),
  OTHER_ALLELE = as.character(native_loadings[["NONMAJ"]]),
  stringsAsFactors = FALSE
)
output <- cbind(
  output,
  native_loadings[, paste0("PC", seq_len(n_pcs)), drop = FALSE]
)

expected_snps <- scan(
  expected_snp_file,
  what = character(),
  quiet = TRUE
)
write_pca_loadings(
  output,
  outfile = output_file,
  n_pcs = n_pcs,
  expected_snps = expected_snps
)

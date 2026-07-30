args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop(
    "Usage: validate_pca_loadings_file.R <pca.loadings.tsv.gz> <n_pcs>"
  )
}

input_file <- args[1]
n_pcs <- as.integer(args[2])

if (!file.exists(input_file) || file.info(input_file)$size == 0L) {
  stop("PCA SNP-loading file is absent or empty: ", input_file)
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
  stop("Could not determine the validation script location.")
}
script_file <- sub("^--file=", "", script_argument)
source(file.path(dirname(normalizePath(script_file)), "pca_loadings.R"))

loadings <- utils::read.delim(
  input_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
validate_pca_loadings(loadings, n_pcs)

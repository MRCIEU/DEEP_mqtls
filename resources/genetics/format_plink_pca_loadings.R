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
pc_columns <- paste0("PC", seq_len(n_pcs))
required_columns <- c(
  chromosome_column, "POS", "ID", "REF", "ALT", "MAJ", "NONMAJ",
  pc_columns
)
missing_columns <- setdiff(required_columns, colnames(native_loadings))
if (length(missing_columns) > 0L) {
  stop(
    "Missing columns in PLINK variant-weight file: ",
    paste(missing_columns, collapse = ", ")
  )
}

if (anyDuplicated(native_loadings[["ID"]])) {
  stop(
    "PLINK variant-weight file contains duplicate variant IDs. ",
    "DEEP expects one row per biallelic variant."
  )
}
if (any(grepl(",", native_loadings[["ALT"]], fixed = TRUE)) ||
    any(grepl(",", native_loadings[["NONMAJ"]], fixed = TRUE))) {
  stop(
    "PLINK variant-weight file contains multi-allelic variant values. ",
    "DEEP expects biallelic variants only."
  )
}

position <- native_loadings[["POS"]]
position <- suppressWarnings(as.integer(position))
if (anyNA(position) || any(position < 1L)) {
  stop("POS values must be positive integers.")
}

# PLINK2 note: with the 'biallelic-var-wts' modifier, an old-style .eigenvec.var file is generated. It's a text file with a header line, followed by one line per variant with the following columns: 

# Header	Column set	Contents
# CHROM	chrom	Chromosome code
# POS	pos	Base-pair coordinate
# ID	(required)	Variant ID
# REF	ref	Reference allele
# ALT1	alt1	Alternate allele 1
# ALT	alt	All alternate alleles, comma-separated
# PROVISIONAL_REF?	maybeprovref, provref	Reports whether REF allele is provisional
# MAJ	maj	Major allele
# NONMAJ	nonmaj	All nonmajor alleles, comma separated
# PC1, PC2, ...	(required)	Principal component variant weights; signs are w.r.t. the major allele

# https://www.cog-genomics.org/plink/2.0/formats

ref_allele <- as.character(native_loadings[["REF"]])
alt_allele <- as.character(native_loadings[["ALT"]])
loading_allele <- as.character(native_loadings[["MAJ"]])
other_allele <- as.character(native_loadings[["NONMAJ"]])
allele_mismatch <- !mapply(
  setequal,
  strsplit(paste(loading_allele, other_allele, sep = "/"), "/", fixed = TRUE),
  strsplit(paste(ref_allele, alt_allele, sep = "/"), "/", fixed = TRUE)
)
if (any(allele_mismatch)) {
  bad_index <- which(allele_mismatch)[1L]
  stop(
    "Could not reconcile MAJ/NONMAJ with REF/ALT for SNP ",
    native_loadings[["ID"]][bad_index],
    ": REF=", ref_allele[bad_index],
    ", ALT=", alt_allele[bad_index],
    ", MAJ=", loading_allele[bad_index],
    ", NONMAJ=", other_allele[bad_index],
    "."
  )
}

output <- data.frame(
  CHR = native_loadings[[chromosome_column]],
  POS = position,
  SNP = as.character(native_loadings[["ID"]]),
  LOADING_ALLELE = loading_allele,
  OTHER_ALLELE = other_allele,
  stringsAsFactors = FALSE
)
output <- cbind(
  output,
  native_loadings[, pc_columns, drop = FALSE]
)

loading_matrix <- as.matrix(output[, pc_columns, drop = FALSE])
storage.mode(loading_matrix) <- "double"
pc_l2 <- sqrt(colSums(loading_matrix^2))
if (any(!is.finite(pc_l2)) || any(pc_l2 <= 0)) {
  bad_pcs <- pc_columns[which(!is.finite(pc_l2) | pc_l2 <= 0)]
  stop(
    "PLINK loading columns have invalid L2 norm for: ",
    paste(bad_pcs, collapse = ", ")
  )
}

normalization_tolerance <- 1e-6
if (all(abs(pc_l2 - 1) <= normalization_tolerance)) {
  message(
    "PLINK loading columns already have unit L2 norm; range before write: ",
    paste(format(range(pc_l2), digits = 8), collapse = " to ")
  )
} else {
  normalized_matrix <- sweep(loading_matrix, MARGIN = 2L, STATS = pc_l2, FUN = "/")
  output[, pc_columns] <- as.data.frame(normalized_matrix)
  normalized_l2 <- sqrt(colSums(normalized_matrix^2))
  message(
    "Normalized PLINK loading columns to unit L2 norm; range before: ",
    paste(format(range(pc_l2), digits = 8), collapse = " to "),
    "; range after: ",
    paste(format(range(normalized_l2), digits = 8), collapse = " to ")
  )
}

allele_audit <- data.frame(
  CHR = native_loadings[[chromosome_column]],
  POS = position,
  SNP = as.character(native_loadings[["ID"]]),
  REF = ref_allele,
  ALT = alt_allele,
  MAJ = loading_allele,
  NONMAJ = other_allele,
  MAJ_IS_REF = loading_allele == ref_allele,
  MAJ_IS_ALT = loading_allele == alt_allele,
  stringsAsFactors = FALSE
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

audit_file <- paste0(output_file, ".allele_audit.tsv")
temporary_audit_file <- tempfile(
  pattern = paste0(".", basename(audit_file), "."),
  tmpdir = dirname(output_file)
)
utils::write.table(
  allele_audit,
  file = temporary_audit_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)
if (!file.rename(temporary_audit_file, audit_file)) {
  copied <- file.copy(temporary_audit_file, audit_file, overwrite = TRUE)
  unlink(temporary_audit_file)
  if (!copied) {
    stop("Could not move the PLINK allele audit table to: ", audit_file)
  }
}
message("Saved PLINK allele audit table to ", audit_file)

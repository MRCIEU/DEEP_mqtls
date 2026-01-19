# Original covariates
# remove genetic PCs
# cell counts
# add HCs

arguments <- commandArgs(T)

comb_cov_file <- arguments[1]
fam_file <- arguments[2]
hcs_file <- arguments[3]
qcovar_hc_file <- arguments[4]
scripts_directory <- arguments[5]

source(paste0(scripts_directory,"/resources/datacheck/fn_rm_constant_col.R"))
source(paste0(scripts_directory,"/resources/datacheck/fn_rm_highlycor.R"))

message("Reading files:")
fam <- read.table(fam_file, header = FALSE)
colnames(fam)[1:2] <- c("FID", "IID")

hcs <- read.table(hcs_file, header = TRUE)
# FID IID HC1 HC2 ... HC100

comb_cov <- read.table(comb_cov_file, header = TRUE, colClass=c("Sex_factor"="character", "Slide_factor"="character"))
# FID IID cell_counts... genetic_pc1 ... genetic_pc20, Age_numeric Sex_factor Slide_factor p_smoking_mcigarette
pcs_exclude <- paste0("genetic_pc", 1:20)
message("Removing genetic PCs from combined covariates file")
comb_cov <- comb_cov[, !(names(comb_cov) %in% pcs_exclude)]

merged <- merge(fam[, 1:2], comb_cov, by = "IID", all.x = TRUE)
if (anyNA(merged)) {
  na_counts <- colSums(is.na(merged))
  na_cols <- na_counts[na_counts > 0]
  stop(sprintf("Merged data contains NA values:\n%s",
               paste(names(na_cols), ":", na_cols, collapse = "\n")))
}

merged <- merge(merged, hcs, by = c("FID", "IID"), all.x = TRUE)
if (anyNA(merged)) {
  na_counts <- colSums(is.na(merged))
  na_cols <- na_counts[na_counts > 0]
  stop(sprintf("Merged data contains NA values after adding HCs:\n%s",
               paste(names(na_cols), ":", na_cols, collapse = "\n")))
}

id_cols <- c("FID", "IID")

factor_cols <- grepl("_factor$", names(merged)) & !(names(merged) %in% id_cols)

merged[factor_cols] <- lapply(merged[factor_cols], as.factor)

numeric_cols <- !(names(merged) %in% id_cols) & !factor_cols

merged[numeric_cols] <- lapply(
  merged[numeric_cols],
  function(x) as.numeric(as.character(x))
)

category_cols <- grepl("_factor$", names(merged)) & !(names(merged) %in% c("FID", "IID"))
quant_cols    <- !(names(merged) %in% c("FID", "IID")) & !category_cols

quant_cov <- merged[, c("FID", "IID", names(merged)[quant_cols & !(names(merged) %in% c("FID", "IID"))])]
category_cov <- merged[, c("FID", "IID", names(merged)[category_cols & !(names(merged) %in% c("FID", "IID"))])]

quant_cov <- remove_constant_cols(quant_cov, "quant_cov")
category_cov <- remove_constant_cols(category_cov, "category_cov")

quant_cov <- remove_constant_cols(quant_cov, "quant_cov")
category_cov <- remove_constant_cols(category_cov, "category_cov")

message("Dectecting prefix in quantitative covariates file")
cell_prefixes <- c("salas\\.", "unilife\\.", "zheng\\.", "middleton\\.")
cell_pattern <- paste0("^(", paste0(cell_prefixes, collapse="|"), ")")
cell_idx <- grep(cell_pattern, colnames(quant_cov))

salas_cols <- grep("^salas\\.", colnames(quant_cov), value = TRUE)
unilife_cols <- grep("^unilife\\.", colnames(quant_cov), value = TRUE)
if (length(salas_cols) > 0 && length(unilife_cols) > 0) {
  message(sprintf("Both 'salas.' and 'unilife.' prefixes detected. Removing salas columns: %s",
                  paste(salas_cols, collapse = ", ")))
  quant_cov <- quant_cov[, !(colnames(quant_cov) %in% salas_cols), drop = FALSE]
  cell_idx <- grep(cell_pattern, colnames(quant_cov))
}

cell_cols <- if (exists("cell_idx") && length(cell_idx) > 0) colnames(quant_cov)[cell_idx] else character(0)

if (length(cell_cols) >= 2) {
  message("Filtering highly correlated cell-type columns (n=", length(cell_cols), "): ", paste(cell_cols, collapse = ", "))
  quant_cov <- filter_correlated_cols(quant_cov, cell_cols, thresh = 0.9, method ="pearson")
} else {
  message("Not enough cell columns for correlation filtering (found ", length(cell_cols), ")")
}

# If any cell-type columns include 'Treg' in their name, remove them to avoid using that column
treg_cols <- grep("treg", colnames(quant_cov), value = TRUE, ignore.case = TRUE)
if (length(treg_cols) > 0) {
  message("Removing columns containing 'Treg' (case-insensitive): ", paste(treg_cols, collapse = ", "))
  quant_cov <- quant_cov[, !(colnames(quant_cov) %in% treg_cols), drop = FALSE]
}

# Check that at least one covariate column remains (besides FID and IID)
if (ncol(quant_cov) < 3) {
  stop("After removing constant columns, quant_cov only contains FID and IID. Please check your input data.")
}
if (ncol(category_cov) < 3) {
  message("After removing constant columns, category_cov only contains FID and IID. Please check your input data. No categorical covariates will be used in GCTA.")
}

qc_cols <- setdiff(colnames(quant_cov), c("FID","IID"))
cc_cols <- setdiff(colnames(category_cov), c("FID","IID"))

qc_uniques <- sapply(quant_cov[, qc_cols, drop = FALSE], function(x) length(unique(x)))
cc_uniques <- sapply(category_cov[, cc_cols, drop = FALSE], function(x) length(unique(as.character(x))))

message("Quantitative covariates unique counts:")
print(qc_uniques)
message("Categorical covariates unique counts:")
print(cc_uniques)

write.table(quant_cov, file = qcovar_hc_file, sep = "\t", row.names = FALSE,  quote = FALSE, col.names = FALSE)

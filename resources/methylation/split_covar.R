# Split covariates file
# Original covariates
# Genetic PCs
# Cell counts

arguments <- commandArgs(T)

comb_cov_file <- arguments[1]
fam_file <- arguments[2]
ccovar_file <- arguments[3]
qcovar_file <- arguments[4]
scripts_directory <- arguments[5]
genetic_pc_gwas <- arguments[6]
qcovar_noPC_file <- arguments[7]

source(paste0(scripts_directory,"/resources/datacheck/fn_rm_constant_col.R"))

for (i in seq_along(arguments)) {
  if (is.na(arguments[i]) || arguments[i] == "") {
    stop(paste0("Argument ", i, " is NA or empty! Please check your input arguments."))
  }
}

# comb_cov file should contain all covariates for mqtl analysis
# including, IID, cell counts prediction, genetic 20 PCs, Age_numeric Sex_factor Slide_factor predicted_smoking
message("Reading files:")
fam <- read.table(fam_file, header = FALSE)
colnames(fam)[1:2] <- c("FID", "IID")

comb_cov <- read.table(comb_cov_file, header = TRUE, colClass=c("Sex_factor"="character", "Slide_factor"="character"))

merged <- merge(fam[, 1:2], comb_cov, by = "IID", all.x = TRUE)

force_category <- c("sex", "Sex", "Sex_factor", "slide", "slide_factor", "Slide_factor")

quant_cols <- sapply(merged, is.numeric)
category_cols <- sapply(merged, function(x) is.character(x) || is.factor(x))

category_cols[names(merged) %in% force_category] <- TRUE
quant_cols[names(merged) %in% force_category] <- FALSE

pcs_exclude <- paste0("genetic_pc", 11:20)
quant_cov <- merged[, c("FID", "IID", names(merged)[quant_cols & !(names(merged) %in% c("FID", "IID", pcs_exclude))])]
category_cov <- merged[, c("FID", "IID", names(merged)[category_cols & !(names(merged) %in% c("FID", "IID"))])]

quant_cov <- remove_constant_cols(quant_cov, "quant_cov")
category_cov <- remove_constant_cols(category_cov, "category_cov")

print(head(category_cov))
print(head(quant_cov))

slide_cols <- grep("Slide_factor", colnames(category_cov), value = TRUE)
if (length(slide_cols) > 0) {
  message("Removing slide factor columns from category_cov: ", paste(slide_cols, collapse = ", "))
  category_cov <- category_cov[, !(colnames(category_cov) %in% slide_cols), drop = FALSE]
}

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

# if (length(cell_idx) > 1) {
#   cell_cols <- colnames(quant_cov)[cell_idx]
#   remove_indices <- integer(0)
#   for (i in seq_along(cell_cols)) {
#     for (j in seq((i+1), length(cell_cols))) {
#       # skip if already marked for removal
#       if (j %in% remove_indices) next
#       x <- as.numeric(quant_cov[[cell_cols[i]]])
#       y <- as.numeric(quant_cov[[cell_cols[j]]])
#       corr_val <- suppressWarnings(cor(x, y, use = "pairwise.complete.obs"))
#       if (!is.na(corr_val) && abs(corr_val) >= 0.9) {
#         # mark the second column (j) for removal
#         remove_indices <- unique(c(remove_indices, j))
#         message(sprintf("Removing collinear cell-count column '%s' (|cor|=%.3f) due to collinearity with '%s'",
#                         cell_cols[j], corr_val, cell_cols[i]))
#       }
#     }
#   }
#   if (length(remove_indices) > 0) {
#     remove_names <- cell_cols[remove_indices]
#     quant_cov <- quant_cov[, !(colnames(quant_cov) %in% remove_names), drop = FALSE]
#   }
# }

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

pc_cols <- paste0("genetic_pc", 1:20)
pc_cov <- quant_cov[, c("FID", "IID", intersect(pc_cols, colnames(quant_cov)))]

all_pc_cols <- paste0("genetic_pc", 1:20)
quant_cov_noPC <- quant_cov[, !(colnames(quant_cov) %in% all_pc_cols), drop = FALSE]

write.table(quant_cov, file = qcovar_file, sep = "\t", row.names = FALSE,  quote = FALSE)
write.table(category_cov, file = ccovar_file, sep = "\t", row.names = FALSE,  quote = FALSE)

write.table(pc_cov, file = genetic_pc_gwas, sep = "\t", row.names = FALSE,  quote = FALSE)
write.table(quant_cov_noPC, file = qcovar_noPC_file, sep = "\t", row.names = FALSE,  quote = FALSE)
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

fam <- read.table(fam_file, header = FALSE)
colnames(fam)[1:2] <- c("FID", "IID")

comb_cov <- read.table(comb_cov_file, header = TRUE, colClass=c("Sex_factor"="character"))

merged <- merge(fam[, 1:2], comb_cov, by = "IID", all.x = TRUE)

force_category <- c("sex", "Sex", "Sex_factor")

quant_cols <- sapply(merged, is.numeric)
category_cols <- sapply(merged, function(x) is.character(x) || is.factor(x))

category_cols[names(merged) %in% force_category] <- TRUE
quant_cols[names(merged) %in% force_category] <- FALSE

quant_cov <- merged[, c("FID", "IID", names(merged)[quant_cols & !(names(merged) %in% c("FID", "IID"))])]
category_cov <- merged[, c("FID", "IID", names(merged)[category_cols & !(names(merged) %in% c("FID", "IID"))])]

quant_cov <- remove_constant_cols(quant_cov, "quant_cov")
category_cov <- remove_constant_cols(category_cov, "category_cov")

# Check that at least one covariate column remains (besides FID and IID)
if (ncol(quant_cov) < 3) {
  stop("After removing constant columns, quant_cov only contains FID and IID. Please check your input data.")
}
if (ncol(category_cov) < 3) {
  stop("After removing constant columns, category_cov only contains FID and IID. Please check your input data.")
}

pc_cols <- paste0("PC", 1:10)
pc_cov <- quant_cov[, c("FID", "IID", intersect(pc_cols, colnames(quant_cov)))]

all_pc_cols <- paste0("PC", 1:20)
quant_cov_noPC <- quant_cov[, !(colnames(quant_cov) %in% all_pc_cols), drop = FALSE]

write.table(quant_cov, file=qcovar_file, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(category_cov, file=ccovar_file, sep = "\t", row.names = FALSE, quote = FALSE)

write.table(pc_cov, file = genetic_pc_gwas, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(quant_cov_noPC, file = qcovar_noPC_file, sep = "\t", row.names = FALSE, quote = FALSE)

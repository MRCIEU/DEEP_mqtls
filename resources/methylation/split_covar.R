# Split covariates file

# Original covariates
# Genetic PCs
# Cell counts

arguments <- commandArgs(T)

comb_cov_file <- arguments[1]
fam_file <- arguments[2]
ccovar_file <- arguments[3]
qcovar_file <- arguments[4]

# comb_cov file should contain all covariates for mqtl analysis
# including, IID, cell counts prediction, genetic 20 PCs, Age_numeric Sex_factor Slide_factor predicted_smoking

fam <- read.table(fam_file, header = FALSE)
colnames(fam)[1:2] <- c("FID", "IID")

comb_cov <- read.table(comb_cov_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

merged <- merge(fam[, 1:2], comb_cov, by = "IID", all.x = TRUE)

force_category <- c("sex", "Sex", "Sex_factor")

quant_cols <- sapply(merged, is.numeric)
category_cols <- sapply(merged, function(x) is.character(x) || is.factor(x))

category_cols[names(merged) %in% force_category] <- TRUE
quant_cols[names(merged) %in% force_category] <- FALSE

quant_cov <- merged[, c("FID", "IID", names(merged)[quant_cols & !(names(merged) %in% c("FID", "IID"))])]
category_cov <- merged[, c("FID", "IID", names(merged)[category_cols & !(names(merged) %in% c("FID", "IID"))])]

write.table(quant_cov, file=qcovar_file, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(category_cov, file=ccovar_file, sep = "\t", row.names = FALSE, quote = FALSE)
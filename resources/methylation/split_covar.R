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

#print(head(category_cov))
#print(head(quant_cov))

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

cell_cols <- if (exists("cell_idx") && length(cell_idx) > 0) colnames(quant_cov)[cell_idx] else character(0)

if (length(cell_cols) >= 2) {
  message("Computing pairwise correlations for cell columns: ", paste(cell_cols, collapse = ", "))
  cor_mat <- tryCatch(
    cor(quant_cov[, cell_cols, drop = FALSE], use = "pairwise.complete.obs"),
    error = function(e) { warning("Correlation matrix failed: ", e$message); return(NULL) }
  )
  if (!is.null(cor_mat)) {
    abs_mat <- abs(cor_mat)
    diag(abs_mat) <- 0
    thresh <- 0.9
    pairs_idx <- which(abs_mat >= thresh, arr.ind = TRUE)
    pairs_idx <- pairs_idx[pairs_idx[,1] < pairs_idx[,2], , drop = FALSE]
    if (nrow(pairs_idx) == 0) {
      message("No cell columns with abs(cor) >= ", thresh)
    } else {
      message("Highly correlated cell column pairs (abs(cor) >= ", thresh, "):")
      apply(pairs_idx, 1, function(ii) {
        i <- ii[1]; j <- ii[2]
        message(sprintf("  %s <-> %s : cor = %.3f", cell_cols[i], cell_cols[j], cor_mat[i,j]))
      })

      # Greedy removal: iteratively remove the column with highest number of high-corr partners;
      # tie-breaker: highest mean absolute correlation with remaining columns.
      to_remove <- character(0)
      remaining <- cell_cols
      adj <- abs_mat >= thresh
      while (TRUE) {
        sub_adj <- adj[remaining, remaining, drop = FALSE]
        if (!any(sub_adj)) break
        deg <- rowSums(sub_adj, na.rm = TRUE)
        maxdeg <- max(deg)
        candidates <- names(deg)[deg == maxdeg]
        if (length(candidates) > 1) {
          meanabs <- rowMeans(abs_mat[remaining, remaining, drop = FALSE], na.rm = TRUE)
          candidate_means <- meanabs[candidates]
          chosen <- candidates[which.max(candidate_means)]
        } else {
          chosen <- candidates[1]
        }
        to_remove <- c(to_remove, chosen)
        message("Removing redundant column: ", chosen, " (degree=", maxdeg, ")")
        remaining <- setdiff(remaining, chosen)
      }

      if (length(to_remove) > 0) {
        quant_cov <- quant_cov[, !(colnames(quant_cov) %in% to_remove), drop = FALSE]
        message("Removed cell columns due to high correlation: ", paste(to_remove, collapse = ", "))
        # update cell_idx to reflect removed columns
        cell_idx <- which(colnames(quant_cov) %in% setdiff(cell_cols, to_remove))
      }
    }
  }
} else {
  message("Not enough cell columns for correlation filtering (found ", length(cell_cols), ")")
}

# If any cell-type columns include 'Treg' in their name, remove them to avoid using that column
treg_cols <- grep("Treg", colnames(quant_cov), value = TRUE)
if (length(treg_cols) > 0) {
  message("Removing columns containing 'Treg': ", paste(treg_cols, collapse = ", "))
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

pc_cols <- paste0("genetic_pc", 1:20)
pc_cov <- quant_cov[, c("FID", "IID", intersect(pc_cols, colnames(quant_cov)))]

all_pc_cols <- paste0("genetic_pc", 1:20)
quant_cov_noPC <- quant_cov[, !(colnames(quant_cov) %in% all_pc_cols), drop = FALSE]

write.table(quant_cov, file = qcovar_file, sep = "\t", row.names = FALSE,  quote = FALSE)
write.table(category_cov, file = ccovar_file, sep = "\t", row.names = FALSE,  quote = FALSE)

write.table(pc_cov, file = genetic_pc_gwas, sep = "\t", row.names = FALSE,  quote = FALSE)
write.table(quant_cov_noPC, file = qcovar_noPC_file, sep = "\t", row.names = FALSE,  quote = FALSE)

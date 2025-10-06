filter_correlated_cols <- function(data, cols, thresh = 0.9, method = c("pearson", "spearman"), verbose = TRUE) {
  method <- match.arg(method)
  if (is.list(cols)) cols <- unlist(cols)
  cols <- as.character(cols)
  cols <- intersect(cols, colnames(data))
  
  if (length(cols) < 2) {
    if (verbose) message("Not enough columns for correlation filtering (found ", length(cols), ")")
    return(data)
  }
  
  cor_mat <- tryCatch(
    cor(data[, cols, drop = FALSE], use = "pairwise.complete.obs", method = method),
    error = function(e) { 
      warning("Correlation matrix failed: ", e$message) 
      return(NULL) 
    }
  )
  if (is.null(cor_mat)) return(data)
  
  abs_mat <- abs(cor_mat)
  diag(abs_mat) <- 0
  pairs_idx <- which(abs_mat >= thresh, arr.ind = TRUE)
  pairs_idx <- pairs_idx[pairs_idx[,1] < pairs_idx[,2], , drop = FALSE]
  if (nrow(pairs_idx) == 0) {
    if (verbose) message("No columns with abs(cor) >= ", thresh)
    return(data)
  }
  
  if (verbose) {
    message("Highly correlated column pairs (abs(cor) >= ", thresh, ", method = ", method, "):")
    apply(pairs_idx, 1, function(ii) {
      i <- ii[1]; j <- ii[2]
      message(sprintf("  %s <-> %s : cor = %.3f", cols[i], cols[j], cor_mat[i,j]))
    })
  }
  
  to_remove <- character(0)
  remaining <- cols
  adj <- abs_mat >= thresh
  while (TRUE) {
    sub_adj <- adj[remaining, remaining, drop = FALSE]
    if (!any(sub_adj, na.rm = TRUE)) break
    deg <- rowSums(sub_adj, na.rm = TRUE)
    maxdeg <- max(deg, na.rm = TRUE)
    candidates <- names(deg)[deg == maxdeg]
    if (length(candidates) > 1) {
      meanabs <- rowMeans(abs_mat[remaining, remaining, drop = FALSE], na.rm = TRUE)
      candidate_means <- meanabs[candidates]
      chosen <- candidates[which.max(candidate_means)]
    } else {
      chosen <- candidates[1]
    }
    to_remove <- c(to_remove, chosen)
    if (verbose) message("Removing redundant column: ", chosen, " (degree=", maxdeg, ")")
    remaining <- setdiff(remaining, chosen)
  }
  
  if (length(to_remove) > 0 && verbose) {
    message("Removed columns due to high correlation: ", paste(to_remove, collapse = ", "))
  }
  
  data[, !(colnames(data) %in% to_remove), drop = FALSE]
}
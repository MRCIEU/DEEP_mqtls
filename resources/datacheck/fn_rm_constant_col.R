remove_constant_cols <- function(df, label, threshold = 0.0003) {
  cols_to_check <- setdiff(names(df), c("FID", "IID"))
  
  col_sds <- sapply(df[cols_to_check], function(x) {
    if (is.numeric(x)) {
      sd(x, na.rm = TRUE)
    } else if (is.factor(x)) {
      if (nlevels(x) <= 1) {
        0
      } else { # leave for T or F
        NA_real_
      }
    } else {
      NA_real_
    }
  })
  
  low_sd_cols <- names(col_sds)[!is.na(col_sds) & col_sds <= threshold]
  
  message("Detected low-variance columns in ", label, ": ",
          ifelse(length(low_sd_cols) > 0, paste(low_sd_cols, collapse = ", "), "none"))
  
  if (length(low_sd_cols) > 0) {
    message(sprintf("In %s, removing columns with SD ≤ %.4f: %s",
                    label, threshold, paste(low_sd_cols, collapse = ", ")))
    df <- df[, !names(df) %in% low_sd_cols, drop = FALSE]
  }
  
  message(sprintf("After filtering, %s columns: %s",
                  label, paste(names(df), collapse = ", ")))
  return(df)
}

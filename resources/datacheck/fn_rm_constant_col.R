remove_constant_cols <- function(df, label) {
  cols_to_check <- setdiff(names(df), c("FID", "IID"))

  constant_cols <- cols_to_check[sapply(df[cols_to_check], function(x) length(unique(x)) == 1)]

  if (length(constant_cols) > 0) {
    message(sprintf("In %s, removing columns with no variation: %s", label, paste(constant_cols, collapse = ", ")))
    df <- df[, !names(df) %in% constant_cols, drop = FALSE]
  }
  
  message(sprintf("After filtering, %s columns: %s", label, paste(names(df), collapse = ", ")))
  return(df)
}
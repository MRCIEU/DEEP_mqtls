#Correlation matrix between measured cell counts and predicted cell counts
library(corrplot)
library(ggplot2)
library(ggExtra)

arguments <- commandArgs(T);
cellcounts_cov <- arguments[1];
measured_cellcounts <- arguments[2];
cellcounts_cov_total <- arguments[3];
cor_matrix <- arguments[4];
cor_plot_ori <- arguments[5];
cor_plot_comb <- arguments[6];
study_name <- arguments[7];
print(study_name)
output_dir <- arguments[8];
print(output_dir)

if (measured_cellcounts != "NULL") {
  message("Reading in measured cell counts")
  # already formatted in 01a
  # colnames should be in prefix of m, e.g. m.Bcells, m.Tcells, m.Mono, m.Gran, m.NK, m.Lym, m.Baso, m.Eos, m.Neu, m.Epi, m.Fib
  # units are in percentage
  measured = read.table(measured_cellcounts, header=T)
  
} else {
  measured = data.frame()
}

message("Reading in predicted cell counts")
predicted<-read.table(cellcounts_cov, header=T)

prefix <- c(unique(sub("\\..*", "", colnames(predicted)[grepl("\\.", colnames(predicted))])))

# standardize predicted column names
standardize_predicted_colnames <- function(data) {
  new_names <- colnames(data)
  # zheng: Neutro -> Neu, Eosino -> Eos
  new_names <- sub("^zheng\\.Neutro$", "zheng.Neu", new_names)
  new_names <- sub("^zheng\\.Eosino$", "zheng.Eos", new_names)
  # middleton: large -> Epi
  new_names <- sub("^middleton\\.large$", "middleton.Epi", new_names)
  colnames(data) <- new_names
  return(data)
}
predicted <- standardize_predicted_colnames(predicted)

print(paste("Found prefixes in predicted cell counts:", paste(prefix, collapse = ", ")))

# there are several situations:
# colnames will be B, CD4T, CD8T, Mono, nRBC (only in babies), Gran, NK, aCD4Tnv, aBaso, aCD4Tmem, aBmem, aBnv, aTreg, aCD8Tmem, aCD8Tnv, aEos, aNK, aNeu, aMono with prefix "unilife."

# or colnames will be CD4Tnv, Baso, CD4Tmem, Bmem, Bnv, Treg, CD8Tmem, CD8Tnv, Eos, NK, Neu, Mono with prefix "salas."
# colnames will be [Epi, Fib, B, NK, CD4T, CD8T, Mono, Neutro] with prefix "zheng."
# colnames will be [CD45pos, large] with prefix "middleton."

# Function to combine cell types for each prefix
combine_cell_types_by_prefix <- function(data, prefixes) {
  combined_data <- data.frame(IID = data$IID)

  for (pref in prefixes) {
    message("Processing prefix: ", pref)
    
    # B cells
    if (pref == "unilife") {
      b_cols <- c(paste0(pref, ".B"), paste0(pref, ".aBmem"), paste0(pref, ".aBnv"))
    } else if (pref == "salas") {
      b_cols <- c(paste0(pref, ".Bmem"), paste0(pref, ".Bnv"))
    } else if (pref == "zheng") {
      b_cols <- c(paste0(pref, ".B"))
    } else {
      b_cols <- character(0)
    }
    available_b_cols <- b_cols[b_cols %in% colnames(data)]
    if (length(available_b_cols) > 0) {
      combined_data[paste0(pref, ".Bcells")] <- rowSums(data[, available_b_cols, drop = FALSE], na.rm = TRUE)
    }

    # T cells
    if (pref == "unilife") {
      t_cols <- c(paste0(pref, ".CD4T"), paste0(pref, ".CD8T"),
                  paste0(pref, ".aCD4Tmem"), paste0(pref, ".aCD4Tnv"), 
                  paste0(pref, ".aTreg"), paste0(pref, ".aCD8Tmem"), paste0(pref, ".aCD8Tnv"))
    } else if (pref == "salas") {
      t_cols <- c(paste0(pref, ".CD4Tmem"), paste0(pref, ".CD4Tnv"), 
                  paste0(pref, ".Treg"), paste0(pref, ".CD8Tmem"), paste0(pref, ".CD8Tnv"))
    } else if (pref == "zheng") {
      t_cols <- c(paste0(pref, ".CD4T"), paste0(pref, ".CD8T"))
    } else {
      t_cols <- character(0)
    }
    available_t_cols <- t_cols[t_cols %in% colnames(data)]
    if (length(available_t_cols) > 0) {
      combined_data[paste0(pref, ".Tcells")] <- rowSums(data[, available_t_cols, drop = FALSE], na.rm = TRUE)
    }
    
    # Monocytes / NK / Gran for unilife
    if (pref == "unilife") {
      mapping <- list(
        "sumMono" = c("Mono", "aMono"),
        "sumNK"   = c("NK", "aNK"),
        "sumGran" = c("Gran", "aNeu", "aEos", "aBaso")
      )
      for (nm in names(mapping)) {
        cols <- paste0(pref, ".", mapping[[nm]])
        available_cols <- cols[cols %in% colnames(data)]
        if (length(available_cols) > 0) {
          combined_data[paste0(pref, ".", nm)] <- rowSums(data[, available_cols, drop = FALSE], na.rm = TRUE)
        }
      }
    }

    # Gran for salas/zheng
    if (pref %in% c("salas", "zheng")) {
      g_cols <- c(paste0(pref, ".Neu"), paste0(pref, ".Eos"), paste0(pref, ".Baso"))
      available_g_cols <- g_cols[g_cols %in% colnames(data)]
      if (length(available_g_cols) > 0) {
        combined_data[paste0(pref, ".Gran")] <- rowSums(data[, available_g_cols, drop = FALSE], na.rm = TRUE)
      }
    }

    if (pref %in% c("unilife","salas", "zheng")) {
      # Lymphocytes = T + B + NK
      t_col <- paste0(pref, ".Tcells")
      b_col <- paste0(pref, ".Bcells")
      nk_col <- if (pref == "unilife") paste0(pref, ".sumNK") else paste0(pref, ".NK")
      available_lymph_cols <- c(t_col, b_col, nk_col)[c(t_col, b_col, nk_col) %in% colnames(combined_data)]
      if (length(available_lymph_cols) > 0) {
        combined_data[paste0(pref, ".Lym")] <- rowSums(combined_data[, available_lymph_cols, drop = FALSE], na.rm = TRUE)
      }
    }
  }
  return(combined_data)
}

# Apply the function to combine cell types by prefix
predicted_comb <- combine_cell_types_by_prefix(predicted, prefix)

message("==== Original predicted columns ====")
message(paste(colnames(predicted), collapse = ", "))

message("==== Final predicted columns after combine_cell_types_by_prefix ====")
message(paste(colnames(predicted_comb), collapse = ", "))

predicted_total <- merge(predicted, predicted_comb, by = "IID", all = TRUE)
write.table(predicted_total, file = cellcounts_cov_total, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

if (nrow(measured) == 0) {
  message("Correlation analysis without measured cell counts")
  # add correlation plot between predicted cell counts
  correlation_matrix1 <- cor(predicted[,-which(names(predicted) == "IID")], use = "complete.obs", method="spearman")

  pdf(file = cor_plot_ori, height = 54, width = 87)
  corrplot(correlation_matrix1, method = "circle", type = "full", tl.col = "black", tl.cex = 5, cl.cex = 5)
  dev.off()

  new_cols <- setdiff(names(predicted_total), names(predicted))
  new_cols <- setdiff(new_cols, "IID")
  correlation_matrix2 <- cor(predicted_total[ , new_cols ], use = "complete.obs", method = "spearman")

  pdf(file = cor_plot_comb, height = 54, width = 87)
  corrplot(correlation_matrix2, method = "circle", type = "full", tl.col = "black", tl.cex = 5, cl.cex = 5)
  dev.off()

  correlation_matrix <- cor(predicted_total[,-which(names(predicted_total) == "IID")], use = "complete.obs", method="spearman")

  data <- predicted_total

} else {
  message("Correlation analysis with measured cell counts")

  data <- merge(measured, predicted_total, by = "IID", all = TRUE)
  ids <- data$IID

  measured <- measured[match(ids, measured$IID), ]
  predicted <- predicted[match(ids, predicted$IID), ]
  
  correlation_matrix <- cor(data[,-which(names(data) == "IID")], use = "complete.obs", method="spearman")

}

write.table(correlation_matrix, file = cor_matrix, row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")

# Scatter plots between predicted cell counts
# Define matching between unilife and salas columns
# unilife: B, CD4T, CD8T, Mono, nRBC (only in babies), Gran, NK from cord blood;
# aCD4Tnv, aBaso, aCD4Tmem, aBmem, aBnv, aTreg, aCD8Tmem, aCD8Tnv, aEos, aNK, aNeu, aMono, sumNK, sumMono, sumEos, sumBaso, sumNeu, sumGran, Bcells, Tcells, Lym
# Tcells, Bcells

# salas: CD4Tnv, Baso, CD4Tmem, Bmem, Bnv, Treg, CD8Tmem, CD8Tnv, Eos, NK, Neu, Mono,
# Tcells, Bcells

# zheng: Epi, Fib, B, NK, CD4T, CD8T, Mono, Neu, Eos, 
# Tcells

# middleton: CD45pos, large

# m: Bcells, Tcells, Lym, Baso, Eos, Neu, Mono, Epi, Fib

method_comp <- list(
  # unilife and salas specific
    CD4Tnv = list(unilife = "unilife.aCD4Tnv", salas = "salas.CD4Tnv"),
    CD4Tmem = list(unilife = "unilife.aCD4Tmem", salas = "salas.CD4Tmem"),
    Bmem = list(unilife = "unilife.aBmem", salas = "salas.Bmem"),
    Bnv = list(unilife = "unilife.aBnv", salas = "salas.Bnv"),
    Treg = list(unilife = "unilife.aTreg", salas = "salas.Treg"),
    CD8Tmem = list(unilife = "unilife.aCD8Tmem", salas = "salas.CD8Tmem"),
    CD8Tnv = list(unilife = "unilife.aCD8Tnv", salas = "salas.CD8Tnv"),

  # unilife, salas, zheng
    B = list(unilife = "unilife.Bcells", salas = "salas.Bcells", m = "m.Bcells", zheng = "zheng.Bcells"),
    T = list(unilife = "unilife.Tcells", salas = "salas.Tcells", m = "m.Tcells", zheng = "zheng.Tcells"),
    Lym = list(unilife = "unilife.Lym", salas = "salas.Lym", m = "m.Lym", zheng = "zheng.Lym"),
    Neu = list(unilife = "unilife.aNeu", salas = "salas.Neu", m= "m.Neu", zheng = "zheng.Neu"),
    Eos = list(unilife = "unilife.aEos", salas = "salas.Eos", m = "m.Eos", zheng = "zheng.Eos"),
    Baso = list(unilife = "unilife.aBaso", salas = "salas.Baso", m = "m.Baso"),
    NK = list(unilife = "unilife.sumNK", salas = "salas.NK", m = "m.NK", zheng = "zheng.NK"),
    Gran = list(unilife = "unilife.sumGran", salas = "salas.Gran", m = "m.Gran"),

    # other cell types, Epi, Fib
    Epi = list(middleton = "middleton.Epi", zheng = "zheng.Epi", m ="m.Epi"),
    Fib = list(zheng = "zheng.Fib", m = "m.Fib")
)

# add correlation (scatter) plots between predicted cell count and between predicted and measured cell counts

compare_and_plot <- function(data, x, y, celltype, prefix_x, prefix_y) {
  if (all(c(x, y) %in% colnames(data))) {
    df <- na.omit(data[, c(x, y)])
    colnames(df) <- c("x", "y")

    if (sd(df$x) == 0 | sd(df$y) == 0) {
      message(paste("Skip plot for", x, "vs", y, ": at least one of them has standard deviation zero"))
      return(NULL)
    }

    corr <- cor(df$x, df$y)
    rmse <- sqrt(mean((df$y - df$x)^2))

    plot_title <- paste0(celltype, ": ", x, " vs ", y)
    suppressMessages(
    p <- ggplot(df, aes(x = x, y = y)) +
      geom_point(alpha = 0.6) +
      geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
      geom_smooth(method = "lm", se = TRUE, color = "blue", linetype = "solid", linewidth = 0.8) +
      labs(
        title = plot_title,
        subtitle = paste0("Pearson r = ", round(corr, 3), ", RMSE = ", round(rmse, 3)),
        x = x,
        y = y
      ) +
      coord_fixed(ratio = 1) + 
      scale_x_continuous(limits = c(0, NA)) +
      scale_y_continuous(limits = c(0, NA)) +
      theme_minimal())

    p_marginal <- ggMarginal(p, type = "density", fill = "lightblue", alpha = 0.5)

    filename <- paste0(output_dir,"/results/01/cellcounts_comp/",study_name,"_scatter_", gsub("\\s+", "_", celltype), "_", prefix_x, "_vs_", prefix_y, ".pdf")
    ggsave(filename, plot = p_marginal, width = 10, height = 9)
  }
}

for (celltype in names(method_comp)) {
  mapping <- method_comp[[celltype]]
  available_vars <- unlist(mapping)
  available_vars <- available_vars[available_vars %in% colnames(data)]
  print(paste("Processing cell type:", celltype, "with variables:", paste(available_vars, collapse = ", ")))

  if (length(available_vars) >= 2) {
    for (i in 1:(length(available_vars)-1)) {
      for (j in (i+1):length(available_vars)) {
        x <- available_vars[i]
        y <- available_vars[j]
        prefix_x <- names(mapping)[which(mapping == x)]
        prefix_y <- names(mapping)[which(mapping == y)]
        compare_and_plot(data, x, y, celltype, prefix_x, prefix_y)
      }
    }
  }
}

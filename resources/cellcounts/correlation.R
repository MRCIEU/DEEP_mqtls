#Correlation matrix between measured cell counts and predicted cell counts
library(corrplot)
library(ggplot2)

arguments <- commandArgs(T);
cellcounts_cov <- arguments[1];
measured_cellcounts <- arguments[2];
cor_matrix <- arguments [3];
cor_plot <- arguments [4];
scripts_directory<-arguments[5];
study_name <- arguments[6];

if (measured_cellcounts != "NULL"){
  message("Reading in measured cell counts")
  # already format in 01a
  # colnames should be in prefix of m, e.g. m.Bcells, m.Tcells, m.Mono, m.Gran, m.NK, m.Lym, m.Baso, m.Eos, m.Neu, m.Epi, m.Fib
  # units are in percentage
  measured = read.table(measured_cellcounts, header=T)
} else {
  measured = data.frame()
}

message("Reading in predicted cell counts")
predicted<-read.table(cellcounts_cov, header=T)
prefix <- unique(sub("\\..*", "", colnames(predicted)[grepl("\\.", colnames(predicted))]))
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
    print(paste("Processing prefix:", pref))
    
    # Get columns for this prefix
    prefix_cols <- colnames(data)[grepl(paste0("^", pref, "\\."), colnames(data))]
    
    if (length(prefix_cols) > 0) {
      print(paste("Using combine_cell_types_by_prefix function, found columns for", pref, ":", paste(prefix_cols, collapse = ", ")))
      
      # B cells combination based on reference type
      if (pref == "unilife") {
        # UniLife: B + aBmem + aBnv
        b_base <- paste0(pref, ".B")
        b_mem <- paste0(pref, ".aBmem")
        b_nv <- paste0(pref, ".aBnv")
        
        available_b_cols <- c(b_base, b_mem, b_nv)[c(b_base, b_mem, b_nv) %in% colnames(data)]
        if (length(available_b_cols) > 0) {
          combined_data[paste0(pref, ".Bcells")] <- rowSums(data[, available_b_cols, drop = FALSE], na.rm = TRUE)
        }
      } else if (pref == "salas") {
        # Salas: Bmem + Bnv
        b_mem <- paste0(pref, ".Bmem")
        b_nv <- paste0(pref, ".Bnv")
        
        available_b_cols <- c(b_mem, b_nv)[c(b_mem, b_nv) %in% colnames(data)]
        if (length(available_b_cols) > 0) {
          combined_data[paste0(pref, ".Bcells")] <- rowSums(data[, available_b_cols, drop = FALSE], na.rm = TRUE)
        }

      } else if (pref == "zheng") {
        # Zheng: B
        b_col <- paste0(pref, ".B")
        if (b_col %in% colnames(data)) {
          combined_data[paste0(pref, ".Bcells")] <- data[, b_col]
        }
      }
      
      # T cells combination based on reference type
      if (pref == "unilife") {
        # UniLife: CD4T + CD8T + aCD4Tmem + aCD4Tnv + aTreg + aCD8Tmem + aCD8Tnv
        t_base <- c(paste0(pref, ".CD4T"), paste0(pref, ".CD8T"))
        t_adult <- c(paste0(pref, ".aCD4Tmem"), paste0(pref, ".aCD4Tnv"), 
                     paste0(pref, ".aTreg"), paste0(pref, ".aCD8Tmem"), paste0(pref, ".aCD8Tnv"))
        
        available_t_cols <- c(t_base, t_adult)[c(t_base, t_adult) %in% colnames(data)]
        if (length(available_t_cols) > 0) {
          combined_data[paste0(pref, ".Tcells")] <- rowSums(data[, available_t_cols, drop = FALSE], na.rm = TRUE)
        }
      } else if (pref == "salas") {
        # Salas: CD4Tmem + CD4Tnv + Treg + CD8Tmem + CD8Tnv
        t_cols <- c(paste0(pref, ".CD4Tmem"), paste0(pref, ".CD4Tnv"), paste0(pref, ".Treg"),
                    paste0(pref, ".CD8Tmem"), paste0(pref, ".CD8Tnv"))
        
        available_t_cols <- t_cols[t_cols %in% colnames(data)]
        if (length(available_t_cols) > 0) {
          combined_data[paste0(pref, ".Tcells")] <- rowSums(data[, available_t_cols, drop = FALSE], na.rm = TRUE)
        }
      } else if (pref == "zheng") {
        # Zheng: CD4T + CD8T
        t_cols <- c(paste0(pref, ".CD4T"), paste0(pref, ".CD8T"))
        
        available_t_cols <- t_cols[t_cols %in% colnames(data)]
        if (length(available_t_cols) > 0) {
          combined_data[paste0(pref, ".Tcells")] <- rowSums(data[, available_t_cols, drop = FALSE], na.rm = TRUE)
        }
      }
      
      # Individual cell types - map to common names
      if (pref == "unilife") {
        # For unilife, combine base and adult cell types
        cell_type_mapping <- list(
          "sumNeu" = c("Neu", "aNeu"),
          "sumMono" = c("Mono", "aMono"),
          "sumEos" = c("Eos", "aEos"),
          "sumBaso" = c("Baso", "aBaso"),
          "sumNK" = c("NK", "aNK"),
          "sumGran" = c("Gran", "aNeu", "aEos", "aBaso") # Granulocytes as sum of Neu, Eos, Baso
        )
        
        for (cell_type in names(cell_type_mapping)) {
          possible_cols <- paste0(pref, ".", cell_type_mapping[[cell_type]])
          available_cols <- possible_cols[possible_cols %in% colnames(data)]
          
          if (length(available_cols) > 0) {
            combined_data[paste0(pref, ".", cell_type)] <- rowSums(data[, available_cols, drop = FALSE], na.rm = TRUE)
          }
        }
      } else {

        # For salas, zheng, middleton - keep original names but standardize Neutro to Neu, large to Epi
        individual_cols <- prefix_cols[!grepl("Bcells|Tcells", prefix_cols)]

        for (col in individual_cols) {
          # Standardize Neutro to Neu for zheng
          if (pref == "zheng" && grepl("Neutro$", col)) {
            new_col_name <- sub("Neutro$", "Neu", col)
            combined_data[new_col_name] <- data[, col]
          } else if (pref == "zheng" && grepl("Eosino$", col)) {
            new_col_name <- sub("Eosino$", "Eos", col)
            combined_data[new_col_name] <- data[, col]
          } else
            combined_data[col] <- data[, col]
          }
        }

        for (col in individual_cols) {
          # standardize large to Epi for middleton
          if (pref == "middleton" && grepl("large$", col)) {
            new_col_name <- sub("large$", "Epi", col)
            combined_data[new_col_name] <- data[, col]
          } else {
            combined_data[col] <- data[, col]
          }
        }

        # Add Granulocytes for salas and zheng
        if (pref %in% c("salas", "zheng")) {
          neu_col <- paste0(pref, ".Neu")
          eos_col <- paste0(pref, ".Eos")
          baso_col <- paste0(pref, ".Baso")
          available_gran_cols <- c(neu_col, eos_col, baso_col)[c(neu_col, eos_col, baso_col) %in% colnames(data)]
          if (length(available_gran_cols) > 0) {
            combined_data[paste0(pref, ".Gran")] <- rowSums(data[, available_gran_cols, drop = FALSE], na.rm = TRUE)
          }
        }
        # no Granulocytes in middleton
      }
    
    # Add Lymphocytes calculation for unilife, salas, zheng (Lymphocytes = Tcells + Bcells + NK)
    for (pref in c("unilife", "salas", "zheng")) {
    t_col <- paste0(pref, ".Tcells")
    b_col <- paste0(pref, ".Bcells")
    # unilife uses sumNK，others use NK
    nk_col <- if (pref == "unilife") paste0(pref, ".sumNK") else paste0(pref, ".NK")
    lymph_col <- paste0(pref, ".Lym")
    available_lymph_cols <- c(t_col, b_col, nk_col)[c(t_col, b_col, nk_col) %in% colnames(predicted)]
    if (length(available_lymph_cols) > 0) {
      predicted[[lymph_col]] <- rowSums(predicted[, available_lymph_cols, drop = FALSE], na.rm = TRUE)
    }
  }
  return(combined_data)
}

# Apply the function to combine cell types by prefix
predicted <- combine_cell_types_by_prefix(predicted, prefix)

print("Final predicted columns:")
print(colnames(predicted))

# Continue with the rest of the correlation analysis

available_cols <- c("IID")
for (pref in prefix) {
  potential_cols <- c(paste0(pref, ".Tcells"), paste0(pref, ".Bcells"), 
                      paste0(pref, ".Neu"), paste0(pref, ".Mono"), 
                      paste0(pref, ".Eos"), paste0(pref, ".Baso"))
  if (pref == "unilife") {
    potential_cols <- c(paste0(pref, ".Tcells"), paste0(pref, ".Bcells"), 
                        paste0(pref, ".sumNeu"), paste0(pref, ".sumMono"), 
                        paste0(pref, ".sumEos"), paste0(pref, ".sumBaso"))
  }
  available_cols <- c(available_cols, potential_cols[potential_cols %in% colnames(predicted)])
}

predicted <- predicted[, colnames(predicted) %in% available_cols]

if (nrow(predicted) == 0) {
  stop("No valid predicted cell counts found.")
}

if (nrow(measured) == 0) {
  print("Correlation analysis without measured cell counts")
  data <- predicted
  # add correlation plot between predicted cell counts
  correlation_matrix <- cor(predicted[,-which(names(predicted) == "IID")], use = "complete.obs", method="spearman")

} else {
  data <- merge(measured, predicted, by = "IID", all = F)
  ids <- data$IID

  measured <- measured[match(ids, measured$IID), ]
  predicted <- predicted[match(ids, predicted$IID), ]
  
  # correlation will include multiple prefixes
  correlation_matrix <- cor(predicted[,-which(names(predicted) == "IID")], measured[,-which(names(measured) == "IID")], use = "complete.obs", method="spearman")

}

write.table(correlation_matrix, file = cor_matrix, row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")
pdf(file = cor_plot, height = 54, width = 87)
corrplot(correlation_matrix, method = "circle", type = "full", tl.col = "black", tl.cex = 5, cl.cex = 5)
dev.off()

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

method_only_comparisons <- list(
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
    Neu = list(unilife = "unilife.sumNeu", salas = "salas.Neu", m= "m.Neu", zheng = "zheng.Neu"),
    Eos = list(unilife = "unilife.sumEos", salas = "salas.Eos", m = "m.Eos", zheng = "zheng.Eos"),
    Baso = list(unilife = "unilife.sumBaso", salas = "salas.Baso", m = "m.Baso"),
    NK = list(unilife = "unilife.sumNK", salas = "salas.NK", m = "m.NK", zheng = "zheng.NK"),
    Gran = list(unilife = "unilife.sumGran", salas = "salas.Gran", m = "m.Gran"),

    # other cell types, Epi, Fib
    Epi = list(middleton = "middleton.Epi", zheng = "zheng.Epi", m ="m.Epi"),
    Fib = list(zheng = "zheng.Fib", m = "m.Fib")
)

# Only plot if both unilife and salas are present
if (all(c("unilife", "salas") %in% prefix)) {
  message("Both unilife and salas predicted cell counts detected. Generating scatter plots between predicted cell counts.")

  for (celltype in names(method_only_comparisons)) {
    unilife_col <- method_only_comparisons[[celltype]]$unilife
    salas_col <- method_only_comparisons[[celltype]]$salas
    if (all(c(unilife_col, salas_col) %in% colnames(predicted))) {
      df <- na.omit(predicted[, c(unilife_col, salas_col)])
      colnames(df) <- c("unilife", "salas")
      corr <- cor(df$unilife, df$salas)
      residuals <- df$salas - df$unilife
      rmse <- sqrt(mean(residuals^2))
      plot_title <- paste0(celltype, ": unilife vs salas")
      p <- ggplot(df, aes(x = unilife, y = salas)) + 
        geom_point(alpha = 0.6) + 
        geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
        geom_smooth(method = "lm", se = TRUE, color = "blue", linetype = "solid", linewidth = 0.8) +
        labs(
          title = plot_title,
          subtitle = paste0("Pearson r = ", round(corr, 3), ", RMSE = ", round(rmse, 3)),
          x = "unilife",
          y = "salas"
        ) +
        coord_fixed(ratio = 1, xlim = c(0, max(c(df$unilife, df$salas), na.rm = TRUE)),
                             ylim = c(0, max(c(df$unilife, df$salas), na.rm = TRUE))) +
        theme_minimal()
      ggsave(paste0(study_name,"_scatter_predicted_", celltype, "_unilife_vs_salas.pdf"), plot = p, width = 6, height = 5)
    }
  }
}

# add correlationn (scatter) plots between predicted cell count and between predicted and measured cell counts


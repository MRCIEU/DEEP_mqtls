#Correlation matrix between measured cell counts and predicted cell counts
library(corrplot)

arguments <- commandArgs(T);
cellcounts_cov <- arguments[1];
measured_cellcounts <- arguments[2];
cor_matrix <- arguments [3];
cor_plot <- arguments [4];
scripts_directory<-arguments[5];

# Define mapping for common cell types and their possible aliases
celltype_aliases <- list(
  Neu = c("neutrophil", "neutrophils", "neutro", "neu", "neut"),
  Lym = c("lymphocyte", "lymphocytes", "lymph", "lym"),
  Mono = c("monocyte", "monocytes", "mono", "mon"),
  Eos = c("eosinophil", "eosinophils", "eos", "eosin"),
  Baso = c("basophil", "basophils", "baso", "bas"),
  Gran = c("granulocyte", "granulocytes", "gran", "grans", "granulo")
)

# Function to standardize measured column names
standardize_measured_colnames <- function(df, mapping, prefix = "m.") {
  cn <- colnames(df)
  for (std in names(mapping)) {
    pats <- mapping[[std]]
    idx <- which(sapply(cn, function(x) any(tolower(x) %in% pats)))
    if (length(idx) == 1) {
      colnames(df)[idx] <- std
      cn <- colnames(df)
    } else if (length(idx) > 1) {
      stop(paste("Multiple columns matched for", std, ":", paste(cn[idx], collapse = ", ")))
    }
  }
  # Add prefix except for IID
  colnames(df)[colnames(df) != "IID"] <- paste0(prefix, colnames(df)[colnames(df) != "IID"])
  return(df)
}

# Function to check if all values (except IID) are <= 1 (assume percent if TRUE)
is_percentage_matrix <- function(df) {
  num_df <- df[, colnames(df) != "IID", drop = FALSE]
  all(num_df <= 1, na.rm = TRUE)
}

# Function to convert to percentage (each row sum to 1)
convert_to_percentage <- function(df) {
  num_df <- df[, colnames(df) != "IID", drop = FALSE]
  row_sums <- rowSums(num_df, na.rm = TRUE)
  # Avoid division by zero
  row_sums[row_sums == 0] <- NA
  num_df <- num_df / row_sums
  df[, colnames(df) != "IID"] <- num_df
  return(df)
}

if (measured_cellcounts != 0) {
  measured <- read.table(measured_cellcounts, header=T)
  if (!"IID" %in% colnames(measured)) {
    stop("Set the name 'IID' to the column of individuals identifiers in measured cell count file")
  }
  print("Detecting columns in measured cell counts file")
  measured <- standardize_measured_colnames(measured, celltype_aliases, prefix = "m.")

  if (nrow(measured) > 0 && !is_percentage_matrix(measured)) {
    message("Measured cell counts are not in percentage, converting to percentage (row sum = 1).")
    measured <- convert_to_percentage(measured)
  }
  
} else {
  measured <- data.frame()
  print("No measured cell counts available for comparison.")
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
        # For salas, zheng, middleton - keep original names but standardize Neutro to Neu
        individual_cols <- prefix_cols[!grepl("Bcells|Tcells", prefix_cols)]
        for (col in individual_cols) {
          # Standardize Neutro to Neu for zheng
          if (pref == "zheng" && grepl("Neutro$", col)) {
            new_col_name <- sub("Neutro$", "Neu", col)
            combined_data[new_col_name] <- data[, col]
          } else {
            combined_data[col] <- data[, col]
          }
        }
      }
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

# add correlation plot between predicted cell count and measured cell counts

# add correlation plot between predicted cell counts

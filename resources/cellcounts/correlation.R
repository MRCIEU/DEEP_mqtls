#Correlation matrix between measured cell counts and predicted cell counts
library(corrplot)

arguments <- commandArgs(T);
cellcounts_cov <- arguments[1];
measured_cellcounts <- arguments[2];
cor_matrix <- arguments [3];
cor_plot <- arguments [4];
scripts_directory<-arguments[5];

measured<-read.table(measured_cellcounts, header=T)

#Check if IID column exist in measured cell counts file 
if (!"IID" %in% colnames(measured)) {
  stop("Set the name 'IID' to the column of individuals identifiers in measured cell count file")
}

# Add prefix 'm.' to all columns except IID
colnames(measured)[colnames(measured) != "IID"] <- paste0("m.", colnames(measured)[colnames(measured) != "IID"])

predicted<-read.table(cellcounts_cov, header=T)
prefix <- unique(sub("\\..*", "", colnames(predicted)[grepl("\\.", colnames(predicted))]))
print(paste("Found prefixes in predicted cell counts:", paste(prefix, collapse = ", ")))

# there are several situations:
# colnames will be B, CD4T, CD8T, Mono, nRBC, Gran, NK, aCD4Tnv, aBaso, aCD4Tmem, aBmem, aBnv, aTreg, aCD8Tmem, aCD8Tnv, aEos, aNK, aNeu, aMono with prefix "unilife."
# QUESTOIN: do we calculate Gran?

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
        # Zheng: already has B
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
          "sumNK" = c("NK", "aNK")
          # add gran if needed
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
data <- merge(measured, predicted, by = "IID", all = F)
ids <- data$IID

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

data <- merge(measured, predicted, by = "IID", all = F)
ids <- data$IID

measured <- measured[match(ids, measured$IID), ]
predicted <- predicted[match(ids, predicted$IID), ]

# correlation will include multiple prefixes
correlation_matrix <- cor(predicted[,-which(names(predicted) == "IID")],
                          measured[,-which(names(measured) == "IID")],
                         use = "complete.obs", method="spearman")

write.table(correlation_matrix, file=cor_matrix, row=TRUE, col=TRUE, qu=FALSE, sep="\t")
pdf(height=54,width=87,cor_plot)
#png(height=540, width=870,file=cor_plot)
corrplot(correlation_matrix, method = "circle", type = "full", tl.col = "black",tl.cex=5,cl.cex=5)
dev.off()
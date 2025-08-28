# Create cell counts file
suppressMessages(library(remotes))
suppressMessages(library(EpiDISH))
suppressMessages(library(meffil))
suppressMessages(library(data.table))

arguments <- commandArgs(T);

tissue<-tolower(arguments[1]);
age<-arguments[2];
methylation_array<-tolower(arguments[3]);

methylation_file<-arguments[4];
cellcounts_cov<-arguments[5];
cellcounts_plot<-arguments[6];
cellcounts_summary<-arguments[7];
scripts_directory<-arguments[8];

# Predict cell counts.
data(cent12CT.m);
data(centUniLIFE.m);
data(centEpiFibIC.m);
data(cent12CT450k.m);
data(centBloodSub.m);

print("Loading methylation file")
load(methylation_file)
print("Sourcing reference selection function")
source(paste0(scripts_directory,"/resources/cellcounts/fn-select_ref.R"))

# Validate inputs and select reference matrices
refs <- validate_and_select_reference(tissue, methylation_array, age)
print(paste0("Using reference matrices ",paste(names(refs), collapse=", ")," for tissue: ", tissue, ", methylation array: ", methylation_array, ", age range: ", age))

cellcounts_total <- data.frame()

for(ref in names(refs)) {
  if (ref %in% c("salas", "unilife")) {
    print(paste0("Using reference: ", ref))
    out.e <- epidish(beta.m = norm.beta, ref.m = refs[[ref]], method = "RPC", maxit = 500)
    cellcounts <- as.data.frame(out.e$estF)
    cellcounts <- as.data.frame(setDT(cellcounts, keep.rownames = "IID"))
    colnames(cellcounts) <- ifelse(
      colnames(cellcounts) != "IID",
       paste0(ref, ".", colnames(cellcounts)),
      colnames(cellcounts)
    )
    # colnames will be B, CD4T, CD8T, Mono, nRBC, Gran, NK, aCD4Tnv, aBaso, aCD4Tmem, aBmem, aBnv, aTreg, aCD8Tmem, aCD8Tnv, aEos, aNK, aNeu, aMono with prefix "unilife."
    # or colnames will be CD4Tnv, Baso, CD4Tmem, Bmem, Bnv, Treg, CD8Tmem, CD8Tnv, Eos, NK, Neu, Mono with prefix "salas."
  } else if (ref == "zheng") {
    print(paste0("Using reference: ", ref))
    out.e <- hepidish(beta.m = norm.beta, ref1.m = centEpiFibIC.m, ref2.m = centBloodSub.m, h.CT.idx = 3, method = 'RPC')
    cellcounts <- as.data.frame(out.e$estF)
    # rownames of IID become the first column
    cellcounts <- as.data.frame(setDT(cellcounts, keep.rownames = "IID"))
    colnames(cellcounts) <- ifelse(
      colnames(cellcounts) != "IID",
       paste0(ref, ".", colnames(cellcounts)),
      colnames(cellcounts)
    )
    # colnames will be [Epi, Fib, B, NK, CD4T, CD8T, Mono, Neutro] with prefix "zheng."
  } else if (ref == "meffil") {
    print(paste0("Using reference: ", ref))
    out.e <- meffil.estimate.cell.counts.from.betas(norm.beta, cell.type.reference = "saliva gse147318")
    cellcounts <- as.data.frame(out.e)
    # rownames of IID become the first column
    cellcounts <- as.data.frame(setDT(cellcounts, keep.rownames = "IID"))
    colnames(cellcounts) <- ifelse(
      colnames(cellcounts) != "IID",
       paste0("middleton.", colnames(cellcounts)),
      colnames(cellcounts)
    )
    # colnames will be [CD45pos, large] with prefix "middleton."
  }

  # Combine results
  if (nrow(cellcounts_total) == 0) {
    cellcounts_total <- cellcounts
  } else {
    cellcounts_total <- merge(cellcounts_total, cellcounts, by = "IID", all = TRUE)
  }
}

#Check if the cellcounts could be calculated for every sample.
cc_na <- table(apply(cellcounts_total, 1, anyNA))
message(sprintf("Cell counts were successfully predicted for %s individuals.", cc_na[1]))
if(length(cc_na) > 1){
  message(sprintf("Cell counts contain NAs for for %s individuals. Please check the input data.", cc_na[2]))
  }

# Loop through each cell type and generate the distribution plot.
pdf(cellcounts_plot, width=12, height=8)
par(mfrow = c(2,3)) # 2x3 grid

for (i in 2:ncol(cellcounts_total)) {
  # Get the cell type for the current iteration
  cell_type <- colnames(cellcounts_total)[i]
  # Generate plots
  plot(cellcounts_total[,i], main = cell_type, xlab = "Sample", ylab = "Estimated Proportion") # Cell count per sample.
  hist(cellcounts_total[,i], main = cell_type, xlab = "Estimated Proportion", ylab = "Frequency", # Histogram of cell counts. 
          col = "lightgrey", border = "black")
  qqnorm(cellcounts_total[,i], xlab = "Theoretical Quantiles", # Predicted cell counts vs. expected distribution of the predicted cell counts.
         ylab = "Estimated Cell Counts", main = cell_type)
  qqline(cellcounts_total[,i], col = "steelblue", lwd = 1)
}

suppressMessages(dev.off())

# Save the distribution of cell counts for this cohort.
library(matrixStats)
cc_mat <- as.matrix(cellcounts_total[-1])
cc_summary <- data.frame(
  mean = colMeans(cc_mat, na.rm = T),
  sd = colSds(cc_mat, na.rm = T),
  min = colMins(cc_mat, na.rm = T),
  perc_0.25 = colQuantiles(cc_mat, probs = 0.25, na.rm = T),
  median = colMedians(cc_mat, na.rm = T),
  perc_0.75 = colQuantiles(cc_mat, probs = 0.75, na.rm = T),
  max = colMaxs(cc_mat, na.rm = T),
  NAs = colSums(is.na(cc_mat))
)

write.table(cc_summary, file = cellcounts_summary, quote = FALSE, row.names = TRUE)

write.table(cellcounts_total, file= cellcounts_cov, row.names=FALSE, col.names=TRUE, quote=FALSE)

message("Cell count prediction complete.\n")
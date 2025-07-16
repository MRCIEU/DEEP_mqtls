# Create cell counts file
library(remotes)
library(EpiDISH)
library(data.table)

arguments <- commandArgs(T)

tissue<- arguments[1]
methylation_file<-arguments[2]
out_file<-arguments[3]
cellcounts_plot<-arguments[4]
cellcounts_summary<-arguments[5]

# Predict cell counts.
data(cent12CT.m)
data(centUniLIFE.m)
load(methylation_file)

out.l1 <- epidish(beta.m = norm.beta, ref.m = cent12CT.m, method = "RPC")
cellcounts<-as.data.frame(out.l1$estF)
cellcounts<-as.data.frame(setDT(cellcounts, keep.rownames = "IID"))
cellcounts <- cellcounts[,c("IID", "Baso", "Bmem", "Bnv", "CD4Tmem", "CD4Tnv", "CD8Tmem", "CD8Tnv", "Eos", "Neu", "NK", "Mono", "Treg")]

out.l1 <- epidish(beta.m = norm.beta, ref.m = centUniLIFE.m, method = "RPC")
cellcounts<-as.data.frame(out.l1$estF)
cellcounts<-as.data.frame(setDT(cellcounts, keep.rownames = "IID"))
cellcounts <- cellcounts.unilife.f[,c("IID", 
                             "B", "CD4T", "CD8T", "Mono", "nRBC", "Gran", "NK", # 7 youthful cord-blood subtypes 
                              "aCD4Tnv", "aBaso", "aCD4Tmem", "aBmem", "aBnv", "aTreg", "aCD8Tmem", 
                              "aCD8Tnv", "aEos", "aNK", "aNeu", "aMono")]

#Check if the cellcounts could be calculated for every sample.
cc_na <- table(apply(cellcounts, 1, anyNA))
message(sprintf("Cell counts were successfully predicted for %s individuals.", cc_na[1]))
if(length(cc_na) > 1){
  message(sprintf("Cell counts contain NAs for for %s individuals. Please check the input data.", cc_na[2]))
  }

# Loop through each cell type and generate the distribution plot.
pdf(cellcounts_plot, width=12, height=8)
par(mfrow = c(2,3)) # 2x3 grid

for (i in 2:ncol(cellcounts)) {
  # Get the cell type for the current iteration
  cell_type <- colnames(cellcounts)[i]
  
  # Generate plots
  plot(cellcounts[,i], main = cell_type, xlab = "Sample", ylab = "Estimated Proportion") # Cell count per sample.
  hist(cellcounts[,i], main = cell_type, xlab = "Estimated Proportion", ylab = "Frequency", # Histogram of cell counts. 
          col = "lightgrey", border = "black")
  qqnorm(cellcounts[,i], xlab = "Theoretical Quantiles", # Predicted cell counts vs. expected distribution of the predicted cell counts.
         ylab = "Estimated Cell Counts", main = cell_type)
  qqline(cellcounts[,i], col = "steelblue", lwd = 1)
}

suppressMessages(dev.off())

# Save the distribution of cell counts for this cohort.
library(matrixStats)
cc_mat <- as.matrix(cellcounts[-1])
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

write.table(cellcounts, file=paste0(out_file), row.names=FALSE, col.names=TRUE, quote=FALSE)

message("Cell count prediction complete.\n")
#!/usr/bin/R

library(GENESIS)

# --- Function: Read GCTA Binary GRM ---
readGRM <- function(rootname) {
  bin.file.name <- paste(rootname, ".grm.bin", sep="")
  n.file.name   <- paste(rootname, ".grm.N.bin", sep="")
  id.file.name  <- paste(rootname, ".grm.id", sep="")
 
  if(!file.exists(bin.file.name)) stop(paste("File not found:", bin.file.name))
  
  id <- read.table(id.file.name)
  n <- dim(id)[1]
  
  bin.file <- file(bin.file.name, "rb")
  grm_val <- readBin(bin.file, n=n*(n+1)/2, what=numeric(0), size=4)
  close(bin.file)
  
  col1 <- rep(1:n, 1:n)
  col2 <- unlist(lapply(1:n, function(i) 1:i))
  
  return(data.frame(id1=col1, id2=col2, grm=grm_val))
}

# --- Command Line Arguments ---
args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 4) stop("Usage: Rscript script.R <gcta_in> <cutoff> <out_prefix> <king_in>")

gcta_infile  <- args[1]
gcta_cutoff  <- as.numeric(args[2])
outfile      <- args[3]
king_infile  <- args[4]

# --- 1. Process GCTA Data ---
message("\n>>> Processing GCTA GRM...")
gcta_data <- readGRM(gcta_infile)
gcta_diag     <- gcta_data$grm[gcta_data$id1 == gcta_data$id2]
gcta_off_diag <- gcta_data$grm[gcta_data$id1 != gcta_data$id2]

# --- 2. Process KING Data ---
message("\n>>> Processing KING Kinship...")
# Ensure king_infile is the correct path to .kin0 or .kin file
king_mat <- kingToMatrix(king = king_infile, estimator = "Kinship")
king_mat <- as.matrix(king_mat)

king_diag     <- diag(king_mat)
king_off_diag <- king_mat[lower.tri(king_mat, diag = FALSE)]
king_cutoff   <- gcta_cutoff / 2

# --- 3. Plotting Top-and-Bottom ---
# Increased height (12) and moderate width (10) for vertical layout
pdf(paste0(outfile, ".pdf"), width = 10, height = 12)
par(mfrow = c(2, 1), mar = c(5, 5, 4, 2)) # Vertical layout with comfortable margins

# Unified Colors and Styles
colors <- c("red", "blue", "cyan", "darkgreen", "purple", "magenta")
ltys <- c(2, 1, 1, 1, 1, 1)

# --- Panel A: GCTA Plot (Top) ---
xmin_g <- -0.2; xmax_g <- 1.2
x_g <- gcta_off_diag
x_plot_g <- x_g[x_g >= xmin_g & x_g <= xmax_g]
h_g <- hist(x_plot_g, breaks = 100, plot = FALSE)
counts_g <- h_g$counts
counts_g[counts_g <= 0] <- 0.5

plot(h_g$mids, counts_g, type = "h", lwd = 3, log = "y", xaxt = "n",
     xlim = c(xmin_g, xmax_g), ylim = c(1, max(counts_g) * 5),
     main = "GCTA GRM Distribution (Off-Diagonal)", 
     xlab = "Genetic Relationship Value", ylab = "Frequency (Log Scale)", col = "darkgray")
axis(1, at = seq(xmin_g, xmax_g, by = 0.1), labels = sprintf("%.1f", seq(xmin_g, xmax_g, by = 0.1)))
abline(v = c(gcta_cutoff, 0, 0.125, 0.25, 0.5, 1.0), col = colors, lwd = 1.5, lty = ltys)

legend("topright", bty = "n", cex = 0.9, title = "GCTA Thresholds",
       legend = c(paste0("Cutoff: ", round(gcta_cutoff, 4)), "0.0: Unrelated", "0.125: 3rd-deg", "0.25: 2nd-deg", "0.5: 1st-deg", "1.0: MZ/Dup"),
       col = colors, lwd = 2, lty = ltys)

# --- Panel B: KING Plot (Bottom) ---
xmin_k <- -0.1; xmax_k <- 0.6
x_k <- king_off_diag
x_plot_k <- x_k[x_k >= xmin_k & x_k <= xmax_k]
h_k <- hist(x_plot_k, breaks = 100, plot = FALSE)
counts_k <- h_k$counts
counts_k[counts_k <= 0] <- 0.5

plot(h_k$mids, counts_k, type = "h", lwd = 3, log = "y", xaxt = "n",
     xlim = c(xmin_k, xmax_k), ylim = c(1, max(counts_k) * 5),
     main = "KING Kinship Distribution (Off-Diagonal)", 
     xlab = "Kinship Coefficient (Phi)", ylab = "Frequency (Log Scale)", col = "darkgray")
axis(1, at = seq(xmin_k, xmax_k, by = 0.05), labels = sprintf("%.2f", seq(xmin_k, xmax_k, by = 0.05)))
abline(v = c(king_cutoff, 0, 0.0625, 0.125, 0.25, 0.5), col = colors, lwd = 1.5, lty = ltys)

legend("topright", bty = "n", cex = 0.9, title = "KING Thresholds",
       legend = c(paste0("Cutoff: ", round(king_cutoff, 4)), "0.0: Unrelated", "0.0625: 3rd-deg", "0.125: 2nd-deg", "0.25: 1st-deg", "0.5: MZ/Dup"),
       col = colors, lwd = 2, lty = ltys)

dev.off()
message("\nDone. Vertical combined plot saved to: ", outfile, ".pdf")
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
# Usage:
#   GCTA only:              Rscript gcta_king_grm_distri.R <gcta_in> <cutoff> <out_prefix>
#   GCTA + KING:            Rscript gcta_king_grm_distri.R <gcta_in> <cutoff> <out_prefix> <kin_in> king
#   GCTA + PC-Relate:       Rscript gcta_king_grm_distri.R <gcta_in> <cutoff> <out_prefix> <kin_in> pcrelate
args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 3) stop("Usage: Rscript gcta_king_grm_distri.R <gcta_in> <cutoff> <out_prefix> [<kin_in> <type>]")

gcta_infile <- args[1]
gcta_cutoff <- as.numeric(args[2])
outfile     <- args[3]
kin_infile  <- if(length(args) >= 4) args[4] else NULL
kin_type    <- if(length(args) >= 5) args[5] else NULL  # "king" or "pcrelate"

has_kin <- !is.null(kin_infile) && !is.null(kin_type)

# --- 1. Process GCTA Data ---
message("\n>>> Processing GCTA GRM...")
gcta_data     <- readGRM(gcta_infile)
message("GRM distribution:")

diag_idx <- gcta_data$grm$id1 == gcta_data$grm$id2
gcta_diagonal <- gcta_data$grm[diag_idx, ]
gcta_off_diag <- gcta_data$grm[!diag_idx, ]

message("\n--- Summary of Diagonal Elements (Self-Relationship) ---")
print(summary(gcta_diagonal$grm))
message("\n--- Summary of Off-Diagonal Elements (Pairwise) ---")
print(summary(gcta_off_diag$grm))

# --- 2. Process Kinship Data (if provided) ---
if (has_kin) {
  if (kin_type == "king") {
    message("\n>>> Processing KING Kinship...")
    king_mat      <- kingToMatrix(king = kin_infile, estimator = "Kinship")
    kin_off_diag  <- as.matrix(king_mat)[lower.tri(as.matrix(king_mat), diag = FALSE)]
    kin_label     <- "KING"
    kin_title     <- "KING Kinship Distribution (Off-Diagonal)"
    kin_xlab      <- "Kinship Coefficient (Phi)"

  } else if (kin_type == "pc_relate") {
    message("\n>>> Processing PC-Relate Kinship...")
    pcrel_data   <- read.table(kin_infile, header = TRUE, sep = "\t")
    kin_off_diag <- pcrel_data$Kinship
    kin_label    <- "PC-Relate"
    kin_title    <- "PC-Relate Kinship Distribution (Pairwise)"
    kin_xlab     <- "Kinship Coefficient (Phi)"

  } else {
    stop(paste("Unknown kin_type:", kin_type, "— must be 'king' or 'pc_relate'"))
  }
  kin_cutoff <- gcta_cutoff / 2
}

# --- 3. Plotting ---
colors <- c("red", "blue", "cyan", "darkgreen", "purple", "magenta")
ltys <- c(2, 1, 1, 1, 1, 1)

if (has_kin) {
  pdf(paste0(outfile, ".pdf"), width = 10, height = 12)
  par(mfrow = c(2, 1), mar = c(5, 5, 4, 2))
} else {
  pdf(paste0(outfile, ".pdf"), width = 10, height = 6)
  par(mfrow = c(1, 1), mar = c(5, 5, 4, 2))
}

# --- Panel A: GCTA ---
xmin_g  <- -0.2; xmax_g <- 1.2
x_plot_g <- gcta_off_diag[gcta_off_diag >= xmin_g & gcta_off_diag <= xmax_g]
h_g     <- hist(x_plot_g, breaks = 100, plot = FALSE)
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

# --- Panel B: KING or PC-Relate (if provided) ---
if (has_kin) {
  xmin_k  <- -0.1; xmax_k <- 0.6
  x_plot_k <- kin_off_diag[kin_off_diag >= xmin_k & kin_off_diag <= xmax_k]
  h_k     <- hist(x_plot_k, breaks = 100, plot = FALSE)
  counts_k <- h_k$counts
  counts_k[counts_k <= 0] <- 0.5

  plot(h_k$mids, counts_k, type = "h", lwd = 3, log = "y", xaxt = "n",
       xlim = c(xmin_k, xmax_k), ylim = c(1, max(counts_k) * 5),
       main = kin_title,
       xlab = kin_xlab, ylab = "Frequency (Log Scale)", col = "darkgray")
  axis(1, at = seq(xmin_k, xmax_k, by = 0.05), labels = sprintf("%.2f", seq(xmin_k, xmax_k, by = 0.05)))
  abline(v = c(kin_cutoff, 0, 0.0625, 0.125, 0.25, 0.5), col = colors, lwd = 1.5, lty = ltys)
  legend("topright", bty = "n", cex = 0.9, title = paste(kin_label, "Thresholds"),
         legend = c(paste0("Cutoff: ", round(kin_cutoff, 4)), "0.0: Unrelated", "0.0625: 3rd-deg", "0.125: 2nd-deg", "0.25: 1st-deg", "0.5: MZ/Dup"),
         col = colors, lwd = 2, lty = ltys)
}

dev.off()
message("\nDone. Plot saved to: ", outfile, ".pdf")

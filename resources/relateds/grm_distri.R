#!/usr/bin/R

#' Read binary GRM files into R
#'
#' @param rootname
#' @export
#' @return List of GRM and id data frames
readGRM <- function(rootname)
{
	bin.file.name <- paste(rootname, ".grm.bin", sep="")
	n.file.name <- paste(rootname, ".grm.N.bin", sep="")
	id.file.name <- paste(rootname, ".grm.id", sep="")
 
	cat("Reading IDs\n")
	id <- read.table(id.file.name)
	n <- dim(id)[1]
	cat("Reading GRM\n")
	bin.file <- file(bin.file.name, "rb")
	grm <- readBin(bin.file, n=n*(n+1)/2, what=numeric(0), size=4)
	close(bin.file)
	cat("Reading N\n")
	n.file <- file(n.file.name, "rb")
	N <- readBin(n.file, n=n*(n+1)/2, what=numeric(0), size=4)
	close(n.file)
 
	cat("Creating data frame\n")
	l <- list()
	for(i in 1:n)
	{
		l[[i]] <- 1:i
	}
	col1 <- rep(1:n, 1:n)
	col2 <- unlist(l)
	grm <- data.frame(id1=col1, id2=col2, N=N, grm=grm)	
 
	ret <- list()
	ret$grm <- grm
	ret$id <- id
	return(ret)
}

arguments <- commandArgs(T)
 
infile <- arguments[1]
cutoff <- as.numeric(arguments[2])
outfile <- arguments[3]

grm <- readGRM(infile)
message("GRM distribution:")

diag_idx <- grm$grm$id1 == grm$grm$id2
grm_diagonal <- grm$grm[diag_idx, ]
grm_off_diag <- grm$grm[!diag_idx, ]

message("\n--- Summary of Diagonal Elements (Self-Relationship) ---")
print(summary(grm_diagonal$grm))
message("\n--- Summary of Off-Diagonal Elements (Pairwise) ---")
print(summary(grm_off_diag$grm))
pdf(paste0(outfile, ".pdf"), width = 10, height = 5)

xmin <- -0.3
xmax <- 1.3
tick <- 0.1
ticks <- seq(xmin, xmax, by = tick)
x <- grm_off_diag$grm

# only keep values within plotting range
x_plot <- x[x >= xmin & x <= xmax]

h <- hist(x_plot, breaks = 100, plot = FALSE)
plot_counts <- h$counts
plot_counts[plot_counts <= 0] <- 0.5 # Use 0.5 to keep it below 1 on log scale

n_lo <- sum(x < xmin, na.rm = TRUE)
n_hi <- sum(x > xmax, na.rm = TRUE)

legend_text <- c(
  paste0("Cutoff = ", cutoff),
  paste0("< ", xmin, ": ", n_lo),
  paste0("> ", xmax, ": ", n_hi),
  "0.00: Unrelated",
  "0.125: Third-degree",
  "0.25: Second-degree",
  "0.50: First-degree",
  "1.00: MZ/Duplicates"
)

colors <- c("red", "black", "black", "blue", "cyan", "darkgreen", "purple", "magenta")

lwds <- c(2, NA, NA, 1, 1, 1, 1, 1)
ltys <- c(2, NA, NA, 1, 1, 1, 1, 1)

plot(
  h$mids, 
  plot_counts, 
  type = "h",             
  lwd = 5,                
  lend = 1,               
  log = "y", 
  xaxt = "n",
  xlim = c(xmin, xmax),
  ylim = c(1, max(plot_counts) * 2),
  main = "GRM Distribution (Log Y-scale)",
  xlab = "Genetic Relationship Value",
  ylab = "Frequency (Log scale)",
  col = "darkgray"
)

axis(1, at = seq(xmin, xmax, by = 0.1), labels = sprintf("%.1f", seq(xmin, xmax, by = 0.1)))

abline(v = cutoff, col = "red", lwd = 2, lty = 2)
abline(v = 0, col = "blue", lwd = 1)
abline(v = 0.125, col = "cyan", lwd = 1)
abline(v = 0.25, col = "darkgreen", lwd = 1)
abline(v = 0.5, col = "purple", lwd = 1)
abline(v = 1.0, col = "magenta", lwd = 1)

legend(
  "topright",
  legend = legend_text,
  col = colors,
  lwd = lwds,
  lty = ltys,
  bty = "n",
  cex = 0.8,
  y.intersp = 1.2
)
dev.off()
message(paste("Plot saved to", outfile, ".pdf"))
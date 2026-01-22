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
print(summary(grm$grm$grm))

pdf(paste0(outfile, ".pdf"), width = 10, height = 5)

xmin <- -0.3
xmax <- 1.3
tick <- 0.1
ticks <- seq(xmin, xmax, by = tick)
x <- grm$grm$grm

# only keep values within plotting range
x_plot <- x[x >= xmin & x <= xmax]

hist(
  x_plot,
  breaks = 500,
  xaxt = "n",
  xlim = c(xmin, xmax),
  main = "GRM distribution",
  xlab = "GRM value"
)

axis(
  1,
  at = ticks,
  labels = sprintf("%.2f", ticks)
)

abline(v = cutoff, col = "red", lwd = 2)

# out-of-range counts
n_lo <- sum(x < xmin, na.rm = TRUE)
n_hi <- sum(x > xmax, na.rm = TRUE)

legend(
  "topright",
  legend = c(
    paste0("cutoff = ", cutoff),
    paste0("< ", xmin, ": ", n_lo),
    paste0("> ", xmax, ": ", n_hi)
  ),
  col = c("red", "black", "black"),
  lwd = c(2, NA, NA),
  bty = "n"
)

dev.off()

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
message(paste0("GRM distribution: ", summary(grm$grm$grm)))

pdf(paste0(outfile, ".pdf"), width = 7, height = 5)
hist(x = grm$grm$grm, breaks = 500,
     main = "GRM distribution",
     xlab = "GRM value")

abline(v = cutoff, col = "red", lwd = 2)
# abline(v = 1, col = "blue", lwd = 2, lty = 2)
legend("topright",
       legend = c(paste0("cutoff=", cutoff)),
       col = c("red"), lwd = 2, lty = c(1), bty = "n")
dev.off()
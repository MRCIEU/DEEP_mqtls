suppressMessages(library(SNPRelate))
suppressMessages(library(GENESIS))
suppressMessages(library(GWASTools))

arguments <- commandArgs(T)
bfile <- arguments[1]
outfile <- arguments[2]
npc <- as.numeric(arguments[3])
nthreads <- as.numeric(arguments[4])
loadings_file <- if(length(arguments) >= 5) {
  arguments[5]
} else {
  paste0(outfile, ".loadings.tsv.gz")
}


message("Reading in genetic data")
bed.fn <- paste0(bfile, ".bed")
bim.fn <- paste0(bfile, ".bim")
fam.fn <- paste0(bfile, ".fam")
gds.fn <- paste0(bfile, ".gds")
snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, gds.fn)


message("Calculating KING-based kinships")
fam <- read.table(fam.fn)
fids <- fam$V1
iids <- fam$V2
genofile <- snpgdsOpen(gds.fn)
ibd <- snpgdsIBDKING(genofile, num.thread=nthreads, family.id=fids)
ibd <- ibd$kinship
colnames(ibd) <- iids
rownames(ibd) <- iids
snpgdsClose(genofile)


message("Calculating PCs based on relatedness")
geno <- GdsGenotypeReader(filename=gds.fn)
genoData <- GenotypeData(geno)
#mypcair <- pcair(genoData = genoData, v=npc, kinMat = ibd, divMat = ibd)
mypcair <- pcair(genoData, kinobj = ibd, divobj = ibd, eigen.cnt = npc)
close(geno)

genofile <- snpgdsOpen(gds.fn)
pca_snp_ids <- gdsfmt::read.gdsn(
  gdsfmt::index.gdsn(genofile, "snp.id")
)
snpgdsClose(genofile)

script_argument <- grep(
  "^--file=",
  commandArgs(trailingOnly = FALSE),
  value = TRUE
)
if (length(script_argument) != 1L) {
  stop("Could not determine the PC-AiR script location.")
}
script_file <- sub("^--file=", "", script_argument)
source(file.path(dirname(normalizePath(script_file)), "pca_loadings.R"))
save_pcair_loadings(
  gds_file = gds.fn,
  mypcair = mypcair,
  pca_snp_ids = pca_snp_ids,
  n_pcs = npc,
  nthreads = nthreads,
  outfile = loadings_file
)

pcs <- mypcair$vectors
ids <- fam[match(fam$V2, rownames(pcs)), 1:2]


message("Saving data")
all(ids$V2 == rownames(pcs))
pcs <- data.frame(ids, pcs)
write.table(pcs, file=paste0(outfile, ".eigenvec"), row=F, col=F, qu=F)
write.table(mypcair$values, file=paste0(outfile, ".eigenval"), row=F, col=F, qu=F)

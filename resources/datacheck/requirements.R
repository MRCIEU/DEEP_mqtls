arguments <- commandArgs(T)
related <- as.character(arguments[1])

message("Checking R version")
currentr <- paste0(R.Version()["major"], ".", R.Version()["minor"])
ch <- compareVersion(currentr, "4.0")
if (ch == -1) {
    stop("You are running R version ", currentr, ". Please upgrade to at least 4.0.")
}

message("Checking that all required packages are present")

pkglist <- c(
    "GwasDataImport",
    "GenomicRanges",
    "lattice",
    "ggplot2",
    "data.table",
    "MatrixEQTL",
    "parallel",
    "matrixStats",
    "plyr",
    "Cairo",
    "plotrix",
    "meffil",
    "EasyQC",
    "impute",
    "EpiDISH",
    "ewaff",
    "MASS",
    "dplyr",
    "magrittr",
    "ggrepel",
    "DunedinPACE",
    "RPMM",
    "vioplot",
    "qqman",
    "ggtext",
    "diptest",
    "mixtools",
    "ggpubr",
    "corrplot",
    "janitor",
    "tidyr",
    "readr",
    "scoreInvHap",
    "EpiDISH"
)

index <- pkglist %in% rownames(installed.packages())
if (any(!index)) {
    stop("Before continuing, the following packages need to be installed:\n", paste(pkglist[!index], collapse = "\n"))
} else {
    message("All required packages installed")
}

pkglist_related <- c("SNPRelate", "GENESIS", "GWASTools")

if (related == "yes") {
    message("Checking that all required packages are present for related samples")

    index <- pkglist_related %in% rownames(installed.packages())

    if (any(!index)) {
        stop("Before continuing, the following packages need to be installed:\n", paste(pkglist_related[!index], collapse = "\n"))
    } else {
        message("All required packages for related samples are installed")
    }
}

# Check for required files
files_to_check <- list(
  # list(dir = "./resources/genetics", pattern = "HRC.r1-1.GRCh37.wgs.mac5.sites.tab.cptid.maf001_recoded.gz"),
  # SNP list
  list(dir = "./resources/genetics", pattern = "topmed.GRCh38.f8wgs.pass.nodup.mac5.maf001.tab.snplist.gz"),
  # chain file for liftover from hg19 to hg38
  #list(dir = "./resources/genetics", pattern = "hg19ToHg38.over.chain"),
  
  list(dir = "./resources/bin/hase/data", pattern = "ref-hrc.ref.gz"),
  list(dir = "./resources/bin/hase/data", pattern = "ref-hrc.ref_info.h5")
)

for (file_info in files_to_check) {
  if (length(list.files(file_info$dir, pattern = file_info$pattern)) == 0) {
    stop(sprintf("Before continuing, you need to download %s from the sftp \n", file_info$pattern))
  }
}

library(meffil)

y <- packageVersion("meffil")
if (y < "1.3.8") {
    stop("Meffil warning: please update to latest version")
}

message("All required packages are installed and required files are downloaded \n\n")

library(EpiDISH)
data(cent12CT.m)
if (exists("cent12CT.m") == F) {
    stop("EpiDISH warning: please install EpiDISH version from github (https://github.com/sjczheng/EpiDISH)")
}

library(DunedinPACE)

x <- packageVersion("DunedinPACE")
if (x < "0.99.0") {
    stop("DunedinPACE warning: please install DunedinPACE version from github (https://github.com/danbelsky/DunedinPACE) but with remotes package")
}

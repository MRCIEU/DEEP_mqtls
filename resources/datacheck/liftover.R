errorlist <- list()
warninglist <- list()

library(data.table)
library(GwasDataImport)

args <- (commandArgs(TRUE));
bim_file <- as.character(args[1]);
ori_build <- as.numeric(args[2]);
miss_liftover <- as.character(args[3])

bim_file <- "/user/work/er20212/DEEP_mqtls/input_data/data_filtered"
ori_build <- 37
miss_liftover <- "/user/work/er20212/DEEP_mqtls/processed_data/miss_liftover.txt"

message("Loading bim file: ", bim_file)
bim <- as.data.frame(fread(paste0(bim_file, ".bim")))

message("Number of SNPs: ", nrow(bim))

message("Checking genome build")
accepted_builds <- c(37, 38)
# determine the genome build
build <- GwasDataImport::determine_build_position(pos=bim[, 4], build = c(37, 38, 36))

if (is.na(build) || !(build %in% accepted_builds)) {
  extra_msg <- if (build == 36) {
    " If you are using build 36, please contact Haotian (haotian.tang@bristol.ac.uk)."
  } else {
    ""
  }
  msg <- sprintf("Please check that the bim file is in build 37 or 38. The bim file does not appear to be in these builds.%s", extra_msg)
  errorlist <- c(errorlist, msg)
  warning("ERROR: ", msg)
}

if (build != ori_build) {
  msg <- sprintf("Please check your config file and bim file. Imputed genome build is GRCh%d, but the bim file is GRCh%d", ori_build, build)
  errorlist <- c(errorlist, msg)
  warning("ERROR: ", msg)
}

if (build == 37) {
    message("Genome build liftover")
    if(mean(grepl("^rs", bim$V3)) > 0.8) {
        temp_bim = GwasDataImport::liftover_gwas(
        dat = bim,
        build = c(37, 38, 36),
        to = 38,
        chr_col = "V1",
        pos_col = "V4",
        snp_col = "V2",
        ea_col = "V5",
        oa_col = "V6")
    } else {
        temp_bim = GwasDataImport::liftover_gwas(
        dat = bim,
        build = c(37, 38, 36),
        to = 38,
        chr_col = "V1",
        pos_col = "V4",
        snp_col = NULL,
        ea_col = "V5",
        oa_col = "V6")
    }

    message(paste0("Liftover complete. ", nrow(temp_bim), " SNPs were successfully converted to GRCh38."))
    message(paste0("Liftover failed for ", nrow(bim) - nrow(temp_bim), " SNPs."))

    if((nrow(temp_bim)/nrow(bim)) < 0.99) {
        msg <- paste0("Liftover was not successful for more than 1% of the SNPs. Please contact Haotian (haotian.tang@bristol.ac.uk).")
        errorlist <- c(errorlist, msg)
        warning("ERROR: ", msg)
    }  
    message("Missing SNPs after liftover saved")
    write.table(bim[!(bim$V2%in%temp_bim$V2),"V2"], file = miss_liftover, sep = "\t", quote = F, row.names = F, col.names = F)

    temp_bim$V2 <- with(temp_bim, paste0(V1, ":", V4, "_",
                        ifelse(V5 < V6, paste(V5, V6, sep = "_"), paste(V6, V5, sep = "_"))))
    bim <- temp_bim
    print(head(temp_bim))
    bim[["V1"]][bim[["V1"]] == "X"] <- "23"
    bim[["V1"]][bim[["V1"]] == "Y"] <- "24"

    write.table(bim, file = paste0(bim_file,"_liftover.bim"), sep = "\t", quote = F, row.names = F, col.names = F)
    } else if (build == 38) {
    message("Genome build is GRCh38, no liftover required")
}

errorlist <- list()
warninglist <- list()

library(data.table)
library(dplyr)
library(GwasDataImport)
library(rtracklayer)

args <- (commandArgs(TRUE));
bim_file <- as.character(args[1]);
ori_build <- as.numeric(args[2]);
miss_liftover <- as.character(args[3])
liftover_map <- as.character(args[4])
section_01_dir <- as.character(args[5])

message("Loading bim file: ", bim_file)
bim <- as.data.frame(fread(paste0(bim_file, ".bim")))

message("Number of SNPs: ", nrow(bim))

message("Checking genome build")
accepted_builds <- c(37, 38)
# determine the genome build
if (mean(startsWith(bim[, 2], "rs")) >= 0.8) {
  message("Detected rsID format, using determine_build function")
  build <- GwasDataImport::determine_build(
    rsid = bim[, 2],
    chr = bim[, 1],
    pos = bim[, 4],
    build = c(37, 38, 36),
    fallback = "position"
  )
} else {
  message("Detected chr:pos format, using determine_build_position function")
  build <- GwasDataImport::determine_build_position(
    pos = bim[, 4],
    build = c(37, 38, 36)
  )
}

write.table(build, file = paste0(section_01_dir, "/inferred_build.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

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
  msg <- sprintf("Please confirm your genome build. Imputed genome build is GRCh%d, but the bim file is GRCh%d", build, ori_build)
  errorlist <- c(errorlist, msg)
  warning("ERROR: ", msg)
}

if (build == 37) {
    message("Genome build liftover")
    if(mean(grepl("^rs", bim$V2)) >= 0.8) {
      # check if the bim file has rsID
      message("Detected rsID format, using liftover_gwas with snp_col")
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
      message("Detected chr:pos format, using liftover_gwas without snp_col")
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
    # temp_bim format: chr, ori_id, V3, new_pos, ref, alt
    message(paste0("Liftover complete. ", nrow(temp_bim), " SNPs were successfully converted to GRCh38."))
    message(paste0("Liftover failed for ", nrow(bim) - nrow(temp_bim), " SNPs."))

    if((nrow(temp_bim)/nrow(bim)) < 0.99) {
        msg <- paste0("Liftover was not successful for more than 1% of the SNPs. Please contact Haotian (haotian.tang@bristol.ac.uk).")
        errorlist <- c(errorlist, msg)
        warning("ERROR: ", msg)
    }

    missing_snps <- bim[!(bim$V2 %in% temp_bim$V2), "V2"]
    if (nrow(as.data.frame(missing_snps)) > 0) {
      message("Save the missing SNPs during liftover")
      write.table(missing_snps, file = miss_liftover, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
      message("Missing SNPs written to file: ", miss_liftover)
      } else {
        message("No missing SNPs during liftover.")
        }

    message("Save the map file")

    map <- left_join(temp_bim, bim, by = "V2")
    map_select <- map[, c(2,4)]
    colnames(map_select) <- c("OLD_SNP", "NEW_POS")

    write.table(map_select, file = liftover_map, sep = "\t", quote = F, row.names = F, col.names = T)

    } else if (build == 38) {

    message("Genome build is GRCh38, no liftover required")

}

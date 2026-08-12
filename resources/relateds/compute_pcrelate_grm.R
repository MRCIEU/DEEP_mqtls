# resources/relateds/compute_pcrelate_grm.R

suppressMessages(library(SNPRelate))
suppressMessages(library(GENESIS))
suppressMessages(library(GWASTools))
suppressMessages(library(Matrix))
if (!requireNamespace("BiocParallel", quietly = TRUE)) BiocManager::install("BiocParallel")
suppressMessages(library(BiocParallel))

args <- commandArgs(trailingOnly = TRUE)
bfile      <- args[1]
outfile    <- args[2]
npc        <- as.numeric(args[3])
nthreads   <- as.numeric(args[4])
rel_cutoff <- as.numeric(args[5])
mode       <- if(length(args) >= 6) args[6] else "keep"  # "keep" or "remove"
loadings_file <- if(length(args) >= 7) {
  args[7]
} else {
  paste0(outfile, ".pca.loadings.tsv.gz")
}
if (length(args) < 8 || !nzchar(args[8]) || args[8] %in% c("NULL", "NA")) {
  stop("HM3 SNP list must be provided as argument 8 for PC-AiR LD pruning.")
}
hm3_snp_list <- args[8]

gds.fn <- paste0(bfile, ".gds")
if (file.exists(gds.fn)) {
  message("Removing existing GDS to avoid stale PC-AiR input: ", gds.fn)
  unlink(gds.fn)
}

snpgdsBED2GDS(paste0(bfile, ".bed"), paste0(bfile, ".fam"), paste0(bfile, ".bim"), gds.fn)

# --- Step 1: SNP Selection (two tracks) ---
genofile <- snpgdsOpen(gds.fn)

# Track A: MAF > 0.01 — for KING kinship + PC-Relate (more SNPs = better precision)
all_snps <- snpgdsSelectSNP(genofile, autosome.only = TRUE, maf = 0.01)

# Track B: MAF > 0.2 + LD pruned — for PC-Air PCA (consistent with plink2 PCA in admix=no)
snps_pca <- snpgdsSelectSNP(genofile, autosome.only = TRUE, maf = 0.2)

if (!file.exists(hm3_snp_list) || file.info(hm3_snp_list)$size == 0L) {
  stop("HM3 SNP list is missing or empty: ", hm3_snp_list)
}
message(">> Restricting PC-AiR LD-pruning candidates to HM3 SNP list: ", hm3_snp_list)
hm3_snp_names <- unique(trimws(readLines(hm3_snp_list, warn = FALSE)))
hm3_snp_names <- hm3_snp_names[nzchar(hm3_snp_names)]
if (length(hm3_snp_names) == 0L) {
  stop("HM3 SNP list contains no SNP IDs: ", hm3_snp_list)
}

gds_snp_id <- read.gdsn(index.gdsn(genofile, "snp.id"))
bim_file <- paste0(bfile, ".bim")
if (!file.exists(bim_file)) {
  stop("Could not find PLINK BIM file for HM3 SNP mapping: ", bim_file)
}
bim <- read.table(bim_file, header = FALSE, stringsAsFactors = FALSE)
if (ncol(bim) < 2L) {
  stop("PLINK BIM file has fewer than two columns: ", bim_file)
}
if (nrow(bim) != length(gds_snp_id)) {
  stop(
    "BIM row count (", nrow(bim), ") does not match GDS snp.id count (",
    length(gds_snp_id), "). Cannot safely map HM3 SNPs."
  )
}
bim_snp_id <- as.character(bim[[2]])
hm3_index <- match(hm3_snp_names, bim_snp_id)
mapped_hm3_snps <- unique(gds_snp_id[hm3_index[!is.na(hm3_index)]])
snps_pca <- intersect(snps_pca, mapped_hm3_snps)
message(
  ">> HM3 list contains ", length(hm3_snp_names), " SNP IDs; ",
  length(mapped_hm3_snps), " map to the GDS; ",
  length(snps_pca), " remain after MAF > 0.2 and autosome filters."
)
if (length(snps_pca) == 0L) {
  stop("No HM3 SNPs remain for PC-AiR LD pruning after mapping and MAF filtering.")
}
plink_r2_threshold <- 0.1
snpgds_ld_threshold <- sqrt(plink_r2_threshold)
message(
  ">> Running SNPRelate LD pruning with method='corr', |r| threshold=",
  signif(snpgds_ld_threshold, 6),
  " to approximate PLINK r^2 threshold=", plink_r2_threshold,
  "."
)
snpset_pruned <- snpgdsLDpruning(
  genofile,
  snp.id = snps_pca,
  method = "corr",
  ld.threshold = snpgds_ld_threshold,
  slide.max.n = 10000,
  num.thread = nthreads,
  start.pos = "first"
)
pruned_snps <- unlist(snpset_pruned)

# --- Step 2: KING with ALL SNPs ---
message(">> Calculating KING matrix with all SNPs...")
king <- snpgdsIBDKING(genofile, num.thread = nthreads, snp.id = all_snps)
king_mat <- king$kinship
colnames(king_mat) <- rownames(king_mat) <- king$sample.id
snpgdsClose(genofile)

# --- Step 3: PC-Air & PC-Relate ---
message(">> Disentangling Ancestry and Kinship...")
geno_reader <- GdsGenotypeReader(filename = gds.fn)
geno_data   <- GenotypeData(geno_reader)

mypcair <- pcair(geno_data, kinobj = king_mat, divobj = king_mat,
                 snp.include = pruned_snps, eigen.cnt = npc)

geno_iter <- GenotypeBlockIterator(geno_data, snpInclude = all_snps)
mypcrel <- pcrelate(geno_iter, pcs = mypcair$vectors[, 1:npc],
                    training.set = mypcair$unrels,
                     BPPARAM = BiocParallel::SerialParam())

# --- Step 4: Export pairwise kinship and GRM ---
kin_out <- data.frame(ID1 = mypcrel$kinBtwn$ID1, ID2 = mypcrel$kinBtwn$ID2, Kinship = mypcrel$kinBtwn$kin)
write.table(kin_out, file = paste0(outfile, ".pcrelate.kin0.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

grm_mat <- pcrelateToMatrix(mypcrel, scaleKin = 2)

grm_sparse <- grm_mat
grm_sparse[grm_sparse < rel_cutoff] <- 0

saveRDS(grm_mat, file = paste0(outfile, ".grm.rds"))
writeMM(as(grm_sparse, "dgCMatrix"), file = paste0(outfile, ".sparse_grm.mtx"))
fam_ids <- read.table(paste0(bfile, ".fam"), header = FALSE, stringsAsFactors = FALSE)
iid <- rownames(mypcair$vectors)
idx <- match(iid, fam_ids$V2)
fid <- fam_ids$V1[idx]
if (any(is.na(fid))) {
  warning("Some IIDs were not found in .fam; using IID as FID for those samples.")
  fid[is.na(fid)] <- iid[is.na(fid)]
}
write.table(data.frame(FID = fid, IID = iid, mypcair$vectors),
            file = paste0(outfile, ".pca.eigenvec"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(mypcair$values,
            file = paste0(outfile, ".pca.eigenval"),
            sep = "\n", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Write GCTA sparse GRM files for fastGWA (--grm-sparse prefix)
message(">> Generating GCTA sparse GRM files for fastGWA...")

# Align grm.id order with matrix row/column order
sample_order <- rownames(grm_sparse)
grm_id_idx <- match(sample_order, iid)
grm_id_df <- data.frame(FID = fid[grm_id_idx], IID = sample_order, stringsAsFactors = FALSE)

write.table(grm_id_df,
            file = paste0(outfile, ".grm.id"),
            sep = " ",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

# pedFAM-style .grm.sp: lower triangle including diagonal, 0-based indices
sp_summary <- summary(as(grm_sparse, "dgCMatrix"))
sp_lower <- sp_summary[sp_summary$i >= sp_summary$j, , drop = FALSE]
sp_lower$i <- sp_lower$i - 1
sp_lower$j <- sp_lower$j - 1

write.table(sp_lower[, c("i", "j", "x")],
            file = paste0(outfile, ".grm.sp"),
            sep = " ",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)
message(">> Written: ", outfile, ".grm.id and ", outfile, ".grm.sp")

# --- Step 5: If mode=remove, identify cryptic relateds to remove ---
if (mode == "remove") {
  message(">> Identifying samples for removal (mode=remove)...")
  kin_threshold <- rel_cutoff / 2  # GRM scale -> kinship scale
  related_pairs <- subset(mypcrel$kinBtwn, kin > kin_threshold)

  remove_iids <- character(0)
  while (nrow(related_pairs) > 0) {
    id_counts <- table(c(related_pairs$ID1, related_pairs$ID2))
    worst     <- names(which.max(id_counts))
    remove_iids <- c(remove_iids, worst)
    related_pairs <- related_pairs[related_pairs$ID1 != worst & related_pairs$ID2 != worst, ]
  }

  fam       <- read.table(paste0(bfile, ".fam"), header = FALSE)
  remove_df <- fam[fam$V2 %in% remove_iids, 1:2]
  write.table(remove_df, file = paste0(outfile, ".remove_ids.txt"),
              sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)
  message(">> ", nrow(remove_df), " samples to remove written to: ", outfile, ".remove_ids.txt")
}

close(geno_reader)

script_argument <- grep(
  "^--file=",
  commandArgs(trailingOnly = FALSE),
  value = TRUE
)
if (length(script_argument) != 1L) {
  stop("Could not determine the PC-Relate script location.")
}
script_file <- sub("^--file=", "", script_argument)
source(file.path(
  dirname(dirname(normalizePath(script_file))),
  "genetics",
  "pca_loadings.R"
))
save_pcair_loadings(
  gds_file = gds.fn,
  mypcair = mypcair,
  pca_snp_ids = pruned_snps,
  n_pcs = npc,
  nthreads = nthreads,
  outfile = loadings_file
)

message(">> PC-Relate Processing finished.")

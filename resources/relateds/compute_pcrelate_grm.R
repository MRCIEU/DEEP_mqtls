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

gds.fn <- paste0(bfile, ".gds")
if (!file.exists(gds.fn)) {
  snpgdsBED2GDS(paste0(bfile, ".bed"), paste0(bfile, ".fam"), paste0(bfile, ".bim"), gds.fn)
}

# --- Step 1: SNP Selection (two tracks) ---
genofile <- snpgdsOpen(gds.fn)

# Track A: MAF > 0.01 — for KING kinship + PC-Relate (more SNPs = better precision)
all_snps <- snpgdsSelectSNP(genofile, autosome.only = TRUE, maf = 0.01)

# Track B: MAF > 0.2 + LD pruned — for PC-Air PCA (consistent with plink2 PCA in admix=no)
snps_pca  <- snpgdsSelectSNP(genofile, autosome.only = TRUE, maf = 0.2)
snpset_pruned <- snpgdsLDpruning(genofile, snp.id = snps_pca, ld.threshold = 0.1,
                                 slide.max.n = 10000, num.thread = nthreads)
pruned_snps <- unlist(snpset_pruned)
write.table(pruned_snps, file = paste0(outfile, ".prune.in"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

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
message(">> PC-Relate Processing finished.")

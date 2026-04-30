# resources/relateds/compute_pcrelate_grm.R

suppressMessages(library(SNPRelate))
suppressMessages(library(GENESIS))
suppressMessages(library(GWASTools))
suppressMessages(library(Matrix))

args <- commandArgs(trailingOnly = TRUE)
bfile      <- args[1]
outfile    <- args[2]
npc        <- as.numeric(args[3])
nthreads   <- as.numeric(args[4])
rel_cutoff <- as.numeric(args[5])

gds.fn <- paste0(bfile, ".gds")
if (!file.exists(gds.fn)) {
  snpgdsBED2GDS(paste0(bfile, ".bed"), paste0(bfile, ".fam"), paste0(bfile, ".bim"), gds.fn)
}

# --- Step 1: SNP Selection ---
genofile <- snpgdsOpen(gds.fn)
# Track A: All SNPs (for high-resolution KING)
all_snps <- snpgdsSelectSNP(genofile, autosome.only = TRUE, maf = 0.01)
# Track B: Pruned SNPs (for robust PCA)
snpset_pruned <- snpgdsLDpruning(genofile, snp.id = all_snps, ld.threshold = 0.2, num.thread = nthreads)
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
                 snps.include = pruned_snps, eigen.cnt = npc)

mypcrel <- pcrelate(geno_data, pcs = mypcair$vectors[, 1:npc], 
                    training.set = mypcair$unrels, snps.include = all_snps)

# --- Step 4: Matrix Processing & Export ---
grm_mat <- pcrelateToMatrix(mypcrel, scaleKin = 2)

grm_sparse <- grm_mat
grm_sparse[grm_sparse < (rel_cutoff * 2)] <- 0

saveRDS(grm_mat, file = paste0(outfile, ".grm.rds"))
writeMM(as(grm_sparse, "dgCMatrix"), file = paste0(outfile, ".sparse_grm.mtx"))
write.table(data.frame(ID = rownames(mypcair$vectors), mypcair$vectors), 
            file = paste0(outfile, ".pca.eigenvec"), sep = "\t", quote = F, row.names = F, col.names = F)

close(geno_reader)
message(">> PC-Relate Processing finished.")
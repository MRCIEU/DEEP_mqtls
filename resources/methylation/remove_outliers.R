library(parallel)
suppressMessages(library(matrixStats))
suppressMessages(library(ewaff))
library(data.table)

arguments <- commandArgs(T);
methylationfile <- arguments[1];
commonids <- arguments[2];
meth_id <- arguments[3];
covar_file<-arguments[4];
fam_file <- arguments[5];
bim_file <- arguments[6];
out_file_gwas <- arguments[7];
out_file_ewas <- arguments[8];
cohort_descriptives_gwas_ewas_file <- arguments[9];
methylation_summary_file_gwas <- arguments[10];
methylation_summary_file_ewas <- arguments[11];
covariates_intersect_gwas <- arguments[12];
covariates_intersect_ewas <- arguments[13];

message("Loading methylation data...")
load(methylationfile)
message("Original methylation input data size: ", ncol(norm.beta), " individuals and ", nrow(norm.beta), " CpGs.")

meth_id<-read.table(meth_id,he=F)
ori_id <- intersect(colnames(norm.beta), meth_id[,1])
norm.beta.meth <- norm.beta[,ori_id]

intersect_ids<-read.table(commonids,he=F)
message("Sample size: ", nrow(intersect_ids), " individuals who have qc-ed genetic data from 01b")

fam <- read.table(fam_file,header=F,stringsAsFactors=F)
shared_samples <- intersect(colnames(norm.beta), fam[,2])
norm.beta.shared <- norm.beta[,shared_samples]
message("Methylation Data size after removal of samples without genetic data: ", ncol(norm.beta.shared), " individuals and ", nrow(norm.beta.shared), " CpGs.")

message("Identifying methylation outliers")

summariseMeth <- function(X, out_file)
{
	message("Counting outliers in methylation matrix")
	
	norm.beta.copy<-ewaff.handle.outliers(X, method=c("iqr"), iqr.limit=3)
	outliers <-norm.beta.copy[[2]]
	keep_cpgs <- rownames(outliers)[which(outliers$n > 0.9*ncol(X))]
	norm.beta <- norm.beta.copy[[1]][rownames(norm.beta.copy[[1]]) %in% keep_cpgs,]
	save(norm.beta, file=out_file)
	message("Data size after removal CpGs with >10% missing values across all samples: ", ncol(norm.beta), " individuals and ", nrow(norm.beta), " CpGs.")
	
	message("Estimating means")
	means <- rowMeans(norm.beta.copy[[1]], na.rm=T)

	message("Estimating SDs")
	sds <- rowSds(norm.beta.copy[[1]], na.rm=T)

	message("Estimating medians")
	medians <- rowMedians(norm.beta.copy[[1]], na.rm=T)
	
	dat <- data.frame(cpg=rownames(norm.beta.copy[[1]]), mean=means, median=medians, sd=sds, outlier=outliers)
	return(dat)
}

message("Generating summary stats of methylation for GWAS")

meth_summary_gwas <- summariseMeth(norm.beta.shared, out_file_gwas)

message("Generating summary stats of methylation for EWAS")

meth_summary_ewas <- summariseMeth(norm.beta.meth, out_file_ewas)

save(meth_summary_gwas, file = methylation_summary_file_gwas)
save(meth_summary_ewas, file = methylation_summary_file_ewas)

message("Reading covariate data and remove outliers")

covs <- read.table(covar_file, header = T, colClasses=c('Sex_factor'='factor'))

message("Data size of covariates: ", nrow(covs), " individuals and ", ncol(covs)-1, " covariates.")

m_gwas<-match(intersect_ids[,1],covs$IID)
covs_gwas<-covs[m_gwas,]
head(covs_gwas)

m_ewas<-match(meth_id[,1],covs$IID)
covs_ewas<-covs[m_ewas,]
head(covs_ewas)

message("Data size of covariates for GWAS: ", nrow(covs_gwas), " individuals and ", ncol(covs_gwas)-1, " covariates.")

message("Data size of covariates for EWAS: ", nrow(covs_ewas), " individuals and ", ncol(covs_ewas)-1, " covariates.")

write.table(covs_gwas,covariates_intersect_gwas,sep="\t",quote=F,row.names=F,col.names=T)
write.table(covs_ewas,covariates_intersect_ewas,sep="\t",quote=F,row.names=F,col.names=T)

bim <- fread(bim_file)

cohort_summary <- list()
cohort_summary$methylation_sample_size_gwas <- ncol(norm.beta.shared)
cohort_summary$methylation_sample_size_ewas <- ncol(norm.beta.meth)
cohort_summary$n_CpGs_gwas <- nrow(norm.beta.shared)
cohort_summary$n_CpGs_ewas <- nrow(norm.beta.meth)
cohort_summary$geno_sample_size <- nrow(fam)
cohort_summary$n_snp <- nrow(bim)
cohort_summary$mqtl_n_males <- sum(covs_gwas$Sex_factor == "M",na.rm=T)
cohort_summary$mqtl_n_females <- sum(covs_gwas$Sex_factor == "F",na.rm=T)
cohort_summary$mqtl_mean_age <- mean(covs_gwas$Age_numeric,na.rm=T)
cohort_summary$mqtl_median_age <- median(covs_gwas$Age_numeric,na.rm=T)
cohort_summary$mqtl_sd_age <- sd(covs_gwas$Age_numeric,na.rm=T)
cohort_summary$mqtl_max_age <- max(covs_gwas$Age_numeric,na.rm=T)
cohort_summary$mqtl_min_age <- min(covs_gwas$Age_numeric,na.rm=T)
cohort_summary$ewas_n_males <- sum(covs_ewas$Sex_factor == "M",na.rm=T)
cohort_summary$ewas_n_females <- sum(covs_ewas$Sex_factor == "F",na.rm=T)
cohort_summary$ewas_mean_age <- mean(covs_ewas$Age_numeric,na.rm=T)
cohort_summary$ewas_median_age <- median(covs_ewas$Age_numeric,na.rm=T)
cohort_summary$ewas_sd_age <- sd(covs_ewas$Age_numeric,na.rm=T)
cohort_summary$ewas_max_age <- max(covs_ewas$Age_numeric,na.rm=T)
cohort_summary$ewas_min_age <- min(covs_ewas$Age_numeric,na.rm=T)
cohort_summary$covariates <- names(covs_gwas)[-1]

save(cohort_summary, file=cohort_descriptives_gwas_ewas_file)

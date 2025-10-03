arguments <- commandArgs(T)

pos_ctr_file <- arguments[1]
filt_pos_ctr_file <- arguments[2]
com_id <- arguments[3]
methylation_no_outliers_gwas <- arguments[4]

library(data.table)
pos_ctr_dt <- fread(pos_ctr_file, header = TRUE)
com_probe_id <- fread(com_id)

message("Check cpg sites from positive control file are all from common probe IDs")

if (all(pos_ctr_dt$positive_control_cpg %in% com_probe_id$Probe_ID)) {
  message("All cpg sites are from common probe IDs")
  filt_pos_ctr_dt <- pos_ctr_dt

} else {

  message(paste0("Remove the following cpg sites are absent from common probe IDs: ", paste(pos_ctr_dt$positive_control_cpg[!pos_ctr_dt$positive_control_cpg %in% com_probe_id$Probe_ID], collapse = ", ")))
  filt_pos_ctr_dt <- pos_ctr_dt[pos_ctr_dt$positive_control_cpg %in% com_probe_id$Probe_ID, ]

}

message("Remove NAs from positive control data")
filt_pos_ctr_dt_nona = na.omit(filt_pos_ctr_dt)
message(paste0("Number of positive control SNP-CpG pairs after filtering NA: ", nrow(filt_pos_ctr_dt_nona)))

load(methylation_no_outliers_gwas)
cpg_names <- rownames(norm.beta)
message("Check cpg sites from positive control file are all in methylation data")
filt_pos_ctr_dt_nona <- filt_pos_ctr_dt_nona[filt_pos_ctr_dt_nona$positive_control_cpg %in% cpg_names, ]
message(paste0("Number of positive control SNP-CpG pairs after filtering CpG sites from meth data: ", nrow(filt_pos_ctr_dt_nona)))

write.table(filt_pos_ctr_dt_nona[,c("positive_control_cpg", "positive_control_snp_chr", "positive_control_snp_pos", "rsid", "positive_control_snp_window", "positive_control_threshold")], file=filt_pos_ctr_file, sep="\t", quote=F, row.names=F, col.names=T)
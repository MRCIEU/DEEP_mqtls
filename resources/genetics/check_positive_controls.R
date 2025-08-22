arguments <- commandArgs(T)

pos_ctr_file <- arguments[1]
filt_pos_ctr_file <- arguments[2]
com_id <- arguments[3]

library(data.table)

pos_ctr_dt <- fread(pos_ctr_file)
com_probe_id <- fread(com_id)

print("Check cpg sites from positive control file are all from common probe IDs")

if (all(pos_ctr_dt$positive_control_cpg %in% com_probe_id$Probe_ID)) {
  print("All cpg sites are from common probe IDs")
  filt_pos_ctr_dt <- pos_ctr_dt
} else {
  print("Remove cpg sites are not from common probe IDs")
  filt_pos_ctr_dt <- pos_ctr_dt[pos_ctr_dt$positive_control_cpg %in% com_probe_id$Probe_ID, ]
}

print("Remove NAs from positive control data")
filt_pos_ctr_dt_nona = na.omit(filt_pos_ctr_dt)
print(paste0("Number of positive control SNP-CpG pairs after filtering: ", nrow(filt_pos_ctr_dt_nona)))

write.table(filt_pos_ctr_dt_nona, file=filt_pos_ctr_file, sep="\t", quote=F, row.names=F, col.names=T)
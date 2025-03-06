library(data.table)

topmed = fread("/user/work/er20212/data/merged_file.snplist", header = TRUE)
hrc = fread("/user/work/er20212/godmc_phase2/resources/genetics/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.cptid.maf001_recoded.gz", header = TRUE)

merged_data = rbind(topmed, hrc)
merged_data$cptid <- sub("^X", "23", merged_data$cptid)
merged_data$cptid <- gsub(" ", "", merged_data$cptid)

merged_data[, c("chr", "pos", "allele1", "allele2") := tstrsplit(cptid, "[:_]")]
merged_data[, al := paste0(allele1, allele2)]
merged_data[, c("allele1", "allele2") := NULL]
merged_data[, `:=`(chr = as.numeric(chr), pos = as.numeric(pos))]
setorder(merged_data, chr, pos, al)
merged_data
undup_merged_data = unique(merged_data, by = "cptid")
undup_merged_data[, c("chr", "pos","al") := NULL]

fwrite(undup_merged_data, "/projects/MRC-IEU/research/projects/ieu3/p3/022/working/data/merged.HRCr1_1.topmedf8.GRCh37.wgs.mac5.maf001.tab.snplist", sep = "\t")

library(data.table)
library(dplyr)

dat = fread("/user/work/er20212/data/pass_filtered.snplist", header = T)

names(dat) <- tolower(names(dat))
head(dat)
nrow(dat)

chunk_size <- 1e6
n <- nrow(dat)
for(i in seq(1, n, by = chunk_size)) {
  end_i <- min(i + chunk_size - 1, n)
  
  dat[i:end_i, ref_alt := ifelse(ref <= alt, paste0(ref, "_", alt),
                                   paste0(alt, "_", ref))]
  dat[i:end_i, cptid := paste0(chr, ":", pos, "_", ref_alt)]
  
  cat("Processed rows", i, "to", end_i, "\n")
}

dat$raf = 1 - dat$af

write.table(dat[,c("cptid","ref","alt","raf")], "/user/work/er20212/data/topmed.GRCh38.f8wgs.pass.mac5.maf001.tab.snplist", sep = "\t", row.names = FALSE, quote = FALSE)

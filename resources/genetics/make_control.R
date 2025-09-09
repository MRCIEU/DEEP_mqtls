library(data.table)
arguments <- commandArgs(T)

in_file <- arguments[1];
ids_file <- arguments[2];
out_file <- arguments[3];

a <- fread(in_file, header=T)
ids <- fread(ids_file, stringsAsFactors=FALSE)
colnames(ids) <- c("FID", "IID")
ids <- ids[match(a$IID, ids$IID), ]
merge <- merge(ids, a, by.x="IID", by.y="IID", sort=F)

write.table(merge, file=out_file, row=F, col=T, qu=F)

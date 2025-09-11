library(data.table)
arguments <- commandArgs(T)

in_file <- arguments[1];
ids_file <- arguments[2];
out_file <- arguments[3];

a <- fread(in_file, header=T)
cpg_id = colnames(a)[2]
ids <- fread(ids_file, stringsAsFactors=FALSE)
colnames(ids) <- c("FID", "IID")
ids <- ids[match(a$IID, ids$IID), ]
merge <- merge(ids, a, by.x="IID", by.y="IID", sort=F)
write.table(merge[, c("FID", "IID", cpg_id), with = FALSE],
            file = out_file,
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)

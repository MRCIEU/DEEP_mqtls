args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]
seed <- as.integer(args[3])

set.seed(seed)
df <- read.csv(input_file, check.names=FALSE)
# 第一列是cpg site，后面是sample id
header <- colnames(df)
df_shuffled <- df
df_shuffled[, -1] <- df[, sample(2:ncol(df), ncol(df)-1)]
write.csv(df_shuffled, output_file, row.names=FALSE, quote=FALSE)
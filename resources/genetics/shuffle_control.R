args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]
seed <- as.integer(args[3])

set.seed(seed)
df <- read.csv(input_file, header = TRUE, sep = "\t")
header <- colnames(df)
df_shuffled <- df
df_shuffled[, 2] <- sample(df[, 2])
write.csv(df_shuffled, output_file, row.names=FALSE, quote=FALSE)
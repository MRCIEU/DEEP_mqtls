args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]
seed <- as.integer(args[3])

set.seed(seed)
df <- read.csv(input_file, header = TRUE, sep = "\t")
header <- colnames(df)
df_shuffled <- df
for (i in 1:nrow(df)) {
  df_shuffled[i, 2:ncol(df)] <- sample(df[i, 2:ncol(df)])
}
write.csv(df_shuffled, output_file, row.names=FALSE, quote=FALSE)
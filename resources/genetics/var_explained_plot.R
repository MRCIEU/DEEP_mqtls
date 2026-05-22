library(ggplot2)
library(data.table)

args        <- commandArgs(trailingOnly = TRUE)
eigenval_file <- args[1]
outfile       <- args[2]

eigenvals <- as.numeric(fread(eigenval_file, header = FALSE)[[1]])
var_exp   <- eigenvals / sum(eigenvals) * 100
df <- data.frame(PC = seq_along(var_exp), var_explained = var_exp)

p <- ggplot(df, aes(x = PC, y = var_explained)) +
    geom_line(colour = "steelblue") +
    geom_point(colour = "steelblue", size = 2) +
    labs(x = "PC", y = "Variance explained (%)") +
    scale_x_continuous(breaks = df$PC) +
    theme_bw() +
    theme(axis.text = element_text(size = 8))

ggsave(plot = p, filename = outfile, height = 4, width = 6)

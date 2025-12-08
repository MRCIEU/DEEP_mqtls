suppressMessages(library(data.table))
library(Matrix)
library(sparsesvd)
suppressMessages(library(ggplot2))

arguments <- commandArgs(T)

full_length_file <- arguments[1]
fam_file <- arguments[2]
output_hcs <- arguments[3]
output_kernel <- arguments[4]
output_plot <- arguments[5]

# no header in fam file
fam <- fread(fam_file, header = FALSE)
sample_id <- data.frame(fam[, .(V1, V2)])
colnames(sample_id) = c("FID", "IID")
nind = nrow(sample_id)

if (nind < 100) {
    number_of_HCs = nind
    message("Number of individuals less than 100, setting number of HCs to ", nind, "\n")
    } else {
    number_of_HCs = 100
    message("Setting number of HCs to 100\n")
}

message("Reading full chunk length data\n")
# no header in the full length file
df <- fread(full_length_file, header = FALSE)

message("Making sparse matrix\n")
A <- sparseMatrix(df$V1, df$V2, x = log10(df$V3+1), dims=c(nind,nind))

message("Performing svd \n")
res<-sparsesvd(A,rank=number_of_HCs)

message("Calculating HCs\n")
HCs <- res$u %*% diag(sqrt(res$d))

message("Writing HCs\n")
colnames(HCs)=paste0("HC",1:number_of_HCs)

# read sample id to add to the first row
HCs <- data.frame(HCs)
HCs <- cbind(sample_id, HCs)
HCs = HCs[, c("FID", "IID", paste0("HC", 1:number_of_HCs))]

fwrite(HCs, output_hcs, row.names = FALSE, sep = ' ', col.names=TRUE)

singvals <- res$d
prop1 <- (singvals) / sum(singvals)
cum_explained  <- cumsum(prop1)
ve <- data.frame(
  HC = paste0("HC", seq_along(singvals)),
  index = seq_along(singvals),
  prop_explained = prop1,
  cum_explained = cum_explained
)

ve_csv <- file.path(output_kernel)
data.table::fwrite(ve, ve_csv)

# Plot: bar of per-HC explained variance, with cumulative line
ve_plot <- ggplot(ve, aes(x = index)) +
  geom_col(aes(y = prop_explained, fill = "prop"), color = NA) +
  geom_line(aes(y = cum_explained, color = "cum prop"), linewidth = 0.9) +
  geom_point(aes(y = cum_explained, color = "cum prop"), size = 1.6) +
  
  # Colors
  scale_fill_manual(values = c("prop" = "#4C78A8"), name = "") +
  scale_color_manual(values = c("cum prop" = "#F58518"), name = "") +
  
  scale_x_continuous(breaks = seq(0, nrow(ve), by = 10)) +
  labs(
    x = "Component index",
    y = "Proportion / Cumulative proportion",
    title = "HC Energy Proportion and Cumulative Proportion"
  ) +
  theme_bw() +
  theme(legend.position = "top")

ve_pdf <- file.path(output_plot)
ggsave(ve_pdf, ve_plot, width = 10, height = 6)
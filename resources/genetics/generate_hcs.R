suppressMessages(library(data.table))
library(Matrix)
library(sparsesvd)
suppressMessages(library(ggplot2))
suppressMessages(library(ggExtra))
suppressMessages(library(gridExtra))

arguments <- commandArgs(T)

full_length_file <- arguments[1]
fam_file <- arguments[2]
output_hcs <- arguments[3]
output_kernel <- arguments[4]
output_plot <- arguments[5]
output_hc_scatter_plot <- arguments[6]

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
ve <- data.frame(
  HC = paste0("HC", seq_along(singvals)),
  index = seq_along(singvals),
  prop_explained = prop1
)

ve_csv <- file.path(output_kernel)
data.table::fwrite(ve, ve_csv)

# Plot: line of per-HC explained variance
ve_plot <- ggplot(ve, aes(x = index, y = prop_explained)) +
  geom_point(size = 1.5) +
  geom_line() +
  scale_y_log10() +
  scale_x_continuous(breaks = ve$index) +
  labs(
    title = "Scree Plot of Top Haplotype Components",
    x     = "HC",
    y     = "Proportion Eigenvalues (log scale)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7)
  )

ve_pdf <- file.path(output_plot)
ggsave(ve_pdf, ve_plot, width = 10, height = 6)

# HC1 vs HCk scatter panel (k = 2..20, capped by available HCs)
# First tile: HC1 marginal density; subsequent tiles: HC1 vs HCk scatter with
# a y-axis marginal density of HCk attached to each subplot.
max_k <- min(20, number_of_HCs)
if (max_k >= 2) {
  hc_df <- as.data.frame(HCs[, c("FID", "IID", paste0("HC", 1:max_k))])

  hc1_density <- ggplot(hc_df, aes(x = HC1)) +
    geom_density(fill = "grey80", colour = "grey40", alpha = 0.6) +
    labs(title = "HC1 distribution", x = "HC1", y = "density") +
    theme_bw() +
    theme(plot.title = element_text(size = 9))

  build_scatter <- function(k) {
    sub <- data.frame(HC1 = hc_df[["HC1"]], HCk = hc_df[[paste0("HC", k)]])
    g <- ggplot(sub, aes(x = HC1, y = HCk)) +
      geom_point(size = 0.6, alpha = 0.6) +
      labs(title = paste0("HC1 vs HC", k), x = "HC1", y = paste0("HC", k)) +
      theme_bw() +
      theme(plot.title = element_text(size = 9))
    ggExtra::ggMarginal(g, margins = "y", type = "density",
               fill = "grey80", colour = "grey40", alpha = 0.6)
  }

  scatter_list <- lapply(2:max_k, build_scatter)
  all_panels <- c(list(hc1_density), scatter_list)

  n_panels <- length(all_panels)
  n_rows   <- ceiling(n_panels / 5)

  combined <- gridExtra::arrangeGrob(grobs = all_panels, ncol = 5,
                          top = "HC1 distribution and HC1 vs HC2..HC20")

  ggsave(output_hc_scatter_plot, combined,
         width = 14, height = max(3 * n_rows, 6), limitsize = FALSE)
} else {
  message("Not enough HCs to plot HC1 vs HCk pairs (need >= 2 HCs)")
}
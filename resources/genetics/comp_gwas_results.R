suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))

args <- commandArgs(T)
filenames <- args[1]
cpg_name <- args[2]
pval_column <- as.numeric(args[3])
beta_column <- as.numeric(args[4])
se_column <- as.numeric(args[5])
chr_column <- as.numeric(args[6])
pos_column  <- as.numeric(args[7])
snp_column  <- as.numeric(args[8])
header <- as.logical(args[9])
control_chr <- as.numeric(args[10])
control_pos <- as.numeric(args[11])
control_window    <- as.numeric(args[12])
control_threshold <- as.numeric(args[13])


filenames <- "/user/work/er20212/alspac_QC_phase1/results/01/01g/positive.control.untransformed.file.txt"
cpg_name <- "cg07959070"
pval_column <- 10
beta_column <- 8
se_column <- 9
chr_column <- 1
pos_column  <- 3
snp_column  <- 2
header <- T
control_chr <- 22
control_pos <- 49660223
control_window    <- 1e5
control_threshold <- 0.001

filenames = read.table(filenames, header = F, sep = "\t")[,1]

gwas_list = list()

for (filename in filenames) {
    message("Reading in ", filename ," GWAS results")
    gwas_result = fread(filename, header = T, data.table=F)
    gwas_result[,chr_column] = as.numeric(gwas_result[,chr_column])
    gwas_result[,pos_column] = as.numeric(gwas_result[,pos_column])
    gwas_result[,pval_column] = as.numeric(gwas_result[,pval_column])
    gwas_result[,snp_column] = as.character(gwas_result[,snp_column])
    gwas_result[,beta_column] = as.numeric(gwas_result[,beta_column])
    gwas_result[,se_column] = as.numeric(gwas_result[,se_column])
    
    index1 <- gwas_result[,chr_column] == control_chr
    index2 <- gwas_result[,pos_column] >= (control_pos - control_window) & gwas_result[,pos_column] <= (control_pos + control_window)
    gwas_result_filter <- gwas_result[index1 & index2, ]
    
    base <- basename(filename)
    outdir <- dirname(filename)

    key <- sub("_.*$", "", base)
    gwas_list[[key]] <- gwas_result_filter
}

keys <- names(gwas_list)
print(keys)
keys <- keys[sapply(gwas_list, nrow) > 0]

if (length(keys) < 2) {
  stop("Less than two non-empty GWAS result sets in gwas_list; cannot compare.")
}

pairs <- combn(keys, 2, simplify = FALSE)

for (pair in pairs) {
  key1 <- pair[1]
  key2 <- pair[2]

  dt1 <- gwas_list[[key1]]
  dt2 <- gwas_list[[key2]]

  df1 <- as.data.frame(dt1)
  df2 <- as.data.frame(dt2)

  snp_column_name = colnames(df1)[snp_column]
  print(snp_column_name)

  merged <- merge(df1[, c(snp_column, beta_column, se_column)],
                df2[, c(snp_column, beta_column, se_column)],
                by = snp_column_name,
                suffixes = c(paste0("_", key1), paste0("_", key2)))

    if (nrow(merged) == 0) {
        message("No overlapping SNPs between ", key1, " and ", key2, "; skip.")
        next
    }

    beta1_name <- paste0(colnames(df1)[beta_column], "_", key1)
    beta2_name <- paste0(colnames(df2)[beta_column], "_", key2)
    se1_name   <- paste0(colnames(df1)[se_column],   "_", key1)
    se2_name   <- paste0(colnames(df2)[se_column],   "_", key2)

  x_col <- beta1_name 
  y_col <- beta2_name
  se_x  <- se1_name    
  se_y  <- se2_name    

  max_abs <- max(abs(merged[[x_col]]), abs(merged[[y_col]]), na.rm = TRUE)
  lim <- c(-max_abs, max_abs)
  brks <- pretty(lim)

  fit <- lm(merged[[y_col]] ~ merged[[x_col]])
  ann_text <- sprintf("slope = %.3f, r = %.3f",
                      coef(fit)[2],
                      cor(merged[[x_col]], merged[[y_col]], use = "complete.obs"))

  p <- ggplot(
    merged,
    aes(x = .data[[x_col]], y = .data[[y_col]])
) +
  geom_errorbarh(
    aes(
      xmin = .data[[x_col]] - .data[[se_x]],
      xmax = .data[[x_col]] + .data[[se_x]]
    ),
    height = 0,
    alpha = 0.4
  ) +
  geom_errorbar(
    aes(
      ymin = .data[[y_col]] - .data[[se_y]],
      ymax = .data[[y_col]] + .data[[se_y]]
    ),
    width = 0,
    alpha = 0.4
  ) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype="dashed", color="red") +
  geom_smooth(method="lm", se=FALSE, color="blue") +
  labs(
    x = paste0("Beta (", key1, ")"),
    y = paste0("Beta (", key2, ")"),
    title    = "Comparison of effect sizes between models",
    subtitle = ann_text
  ) +
  theme_minimal()

  out_pdf <- file.path(outdir, paste0(cpg_name, "_", key1, "_", key2, "_scatter.pdf"))
  ggsave(out_pdf, p, width = 6, height = 6)
  message("Saved scatter plot: ", out_pdf)
}
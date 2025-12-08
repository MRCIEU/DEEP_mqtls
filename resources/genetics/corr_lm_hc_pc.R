## Compute correlations between HCs (hc1..hc100) and PCs (PC1..PC10)
## make scatter plots (one PDF per PC with 10 panels).
## linear model: HC predicts PC (HC -> PC)
## linear model: PC predicts HC (PC -> HC)
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
library(gridExtra)

arguments <- commandArgs(T)
pca_file = arguments[1]
hcs_file = arguments[2]
genetic_outlier_ids = arguments[3]
out_dir = arguments[4]

## Read data
PCs <- fread(pca_file)
names(PCs) <- c("FID","IID", paste("genetic_pc", 1:(ncol(PCs)-2), sep=""))
n_pc_avail <- ncol(PCs) - 2
if (n_pc_avail < 10) {
  stop(sprintf("Expected at least 10 PCs, but found %d. Abort.", n_pc_avail))
}
names(PCs) <- c("FID","IID", paste0("genetic_pc", 1:n_pc_avail))
PCs <- PCs[, c("FID", "IID", paste0("genetic_pc", 1:10))]

# FID, IID, HC1..HCn
HCs = fread(hcs_file)

# remove FID column if present
if("FID" %in% names(HCs)) HCs = HCs[, FID := NULL]
if("FID" %in% names(PCs)) PCs = PCs[, FID := NULL]

merged = merge(HCs, PCs, by = "IID")

hc_cols <- grep("(?i)^hc\\d+$", names(merged), value = TRUE, perl = TRUE)
pc_cols <- grep("^genetic_pc\\d+$", names(merged), value = TRUE)
if (length(hc_cols) == 0) stop("No HC columns found (expected HC1..HCn)")
if (length(pc_cols) == 0) stop("No PC columns found (expected genetic_pc1..genetic_pc10)")
message(sprintf("Found %d HC columns and %d PC columns", length(hc_cols), length(pc_cols)))

## Prepare output containers and directories
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

res = data.table(hc=character(), pc=character(), r=numeric(), p.value=numeric(), n=integer())

for (pc in pc_cols) {
  # Build long data: one row per IID×HC for this PC
  dt <- merged[, c("IID", pc, hc_cols), with = FALSE]
  setnames(dt, pc, "pc_value")
  long <- melt(dt, id.vars = c("IID", "pc_value"),
               measure.vars = hc_cols, variable.name = "HC", value.name = "hc_value")

  # Per-HC stats for this PC
  stats <- long[, {
    ok <- complete.cases(hc_value, pc_value)
    n <- sum(ok)
    if (n < 3) {
      .(r = NA_real_, p.value = NA_real_, n = n, slope = NA_real_)
    } else {
      mod <- lm(pc_value[ok] ~ hc_value[ok])
      ct <- coef(summary(mod))
      r <- cor(pc_value[ok], hc_value[ok], use = "complete.obs")
      pval <- ifelse(nrow(ct) >= 2, ct[2, "Pr(>|t|)"], NA_real_)
      slope <- ifelse(nrow(ct) >= 2, ct[2, "Estimate"], NA_real_)
      .(r = r, p.value = pval, n = n, slope = slope)
    }
  }, by = HC]

  labels_tab <- stats[, .(
  HC,
  lab = sprintf("%s\nr=%.2f  p=%s",
                HC, r, format.pval(p.value, digits = 2))
  )]

  long[, HC_lab := labels_tab$lab[match(HC, labels_tab$HC)]]

  p <- ggplot(long, aes(x = hc_value, y = pc_value)) +
  geom_point(size = 0.6, alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "red", na.rm = TRUE) +
  facet_wrap(~ HC_lab, ncol = 10, scales = "free_y") +
  theme_bw() +
  theme(
    strip.text = element_text(size = 7, lineheight = 0.9),
    strip.background = element_rect(fill = "grey95")
  ) +
  labs(x = "HC value", y = pc, title = sprintf("%s vs all HCs", pc))

  # Dynamic page size based on number of HCs
  n_hc <- length(hc_cols); ncol <- 10; nrow <- ceiling(n_hc / ncol)
  width <- max(12, 1.6 * ncol); height <- max(8, 1.6 * nrow)

  pdf_file <- file.path(out_dir, paste0(pc, ".HCs_scatter.pdf"))
  ggsave(pdf_file, p, width = width, height = height)
  message(sprintf("Wrote %s", pdf_file))

  # Collect stats
  res <- rbind(res, stats[, .(hc = HC, pc = pc, r = r, p.value = p.value, n = n)])
}

# Save correlation results table
fwrite(res, file = file.path(out_dir, "hc_pc_correlations.csv"))


# Linear models: HC predicts PC (HC -> PC) and PC predicts HC (PC -> HC)

get_model_stats <- function(fit, actual, pred) {
  r2  <- summary(fit)$r.squared
  r2_adj <- summary(fit)$adj.r.squared
  
  # Pearson correlation
  r <- cor(actual, pred)
  
  # RMSE
  rmse <- sqrt(mean((actual - pred)^2))
  
  list(r2 = r2, r2_adj = r2_adj, r = r, rmse = rmse)
}

# --------------------------------------------------------------------
# 1) All_PCs_Predicted_vs_Actual.pdf
# --------------------------------------------------------------------
pdf(file.path(out_dir, "All_PCs_Predicted_vs_Actual.pdf"),
    width = 20, height = 10)

plot_list <- list()

for (pc in pc_cols) {

  dt <- merged[, c(pc, hc_cols), with = FALSE]
  setnames(dt, pc, "pc_value")
  dt <- dt[complete.cases(dt)]

  # Model
  formula_str <- paste("pc_value ~", paste(hc_cols, collapse = " + "))
  fit <- lm(as.formula(formula_str), data = dt)

  # Predictions
  dt[, pred := predict(fit)]

  # Stats
  stats <- get_model_stats(fit, dt$pc_value, dt$pred)

  title_str <- sprintf("%s (R² = %.3f, RMSE = %.3f)", pc, stats$r2, stats$rmse)

  p <- ggplot(dt, aes(x = pred, y = pc_value)) +
    geom_point(alpha = 0.6, size = 0.8) +
    geom_smooth(method = "lm", color = "red", se = FALSE) +
    theme_bw(base_size = 12) +
    labs(
      x = "Predicted PC value",
      y = "Actual PC value",
      title = title_str
    )

  plot_list[[pc]] <- p
}

# Arrange plots
n <- length(plot_list)
ncol <- 4
nrow <- ceiling(n / ncol)
do.call("grid.arrange", c(plot_list, ncol = ncol, nrow = nrow))

dev.off()


# --------------------------------------------------------------------
# 2) All_HCs_Predicted_vs_Actual.pdf
# --------------------------------------------------------------------
pdf(file.path(out_dir, "All_HCs_Predicted_vs_Actual.pdf"),
    width = 22, height = 30)  # wide enough for 6 per row

plot_list <- list()

for (hc in hc_cols) {

  dt <- merged[, c(hc, pc_cols), with = FALSE]
  setnames(dt, hc, "hc_value")
  dt <- dt[complete.cases(dt)]

  # Model
  formula_str <- paste("hc_value ~", paste(pc_cols, collapse = " + "))
  fit <- lm(as.formula(formula_str), data = dt)

  # Predictions
  dt[, pred := predict(fit)]

  # Stats
  stats <- get_model_stats(fit, dt$hc_value, dt$pred)

  # Annotation text
  title_str <- sprintf("%s (R² = %.3f, RMSE = %.3f)", hc, stats$r2, stats$rmse)

  p <- ggplot(dt, aes(x = pred, y = hc_value)) +
    geom_point(alpha = 0.6, size = 0.7) +
    geom_smooth(method = "lm", color = "blue", se = FALSE) +
    theme_bw(base_size = 10) +
    labs(
      x = "Predicted HC value",
      y = "Actual HC value",
      title = title_str
    )

  plot_list[[hc]] <- p
}

# 10 plots per row
n <- length(plot_list)
ncol <- 6
nrow <- ceiling(n / ncol)
do.call("grid.arrange", c(plot_list, ncol = ncol, nrow = nrow))

dev.off()

message("Saved both enhanced PC and HC prediction PDFs with all statistics.")

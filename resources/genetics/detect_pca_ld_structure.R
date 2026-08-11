#!/usr/bin/env Rscript

# Detect long-range LD structure from PCA SNP loadings.
#
# NOTE ON METHOD
# --------------
# This script deliberately follows the single-round outlier-detection step used
# by bigsnpr::bed_autoSVD() in Privé et al. (2020):
#
#   1. S <- sqrt(bigutilsr::dist_ogk(loadings))
#   2. Smooth S with bigutilsr::rollmean(), separately by chromosome.
#   3. threshold <- bigutilsr::tukey_mc_up(smoothed_S, alpha = 0.05)
#   4. Flag variants with smoothed_S > threshold.
#   5. Summarise runs of at least 20 consecutive flagged variants as candidate
#      long-range LD regions.
#
# The original bed_autoSVD() procedure removes flagged variants, recomputes the
# PCA, and repeats. DEEP runs this script in 01b as a QC report only, after the
# formal 01b PCs and SNP loadings have been written. It therefore MUST NOT alter
# pca.eigenvec, pca.eigenval, the genotype data, or the loading file.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop(
    "Usage: detect_pca_ld_structure.R ",
    "<pca.loadings.tsv.gz> <output_prefix> <n_report_pcs> ",
    "<roll_radius> <tukey_alpha> <min_region_snps>"
  )
}

loadings_file <- args[1]
output_prefix <- args[2]
n_report_pcs <- as.integer(args[3])
roll_radius <- as.integer(args[4])
tukey_alpha <- as.numeric(args[5])
min_region_snps <- as.integer(args[6])

if (!file.exists(loadings_file) || file.info(loadings_file)$size == 0L) {
  stop("PCA SNP-loading file is absent or empty: ", loadings_file)
}
if (is.na(n_report_pcs) || n_report_pcs < 2L) {
  stop("n_report_pcs must be an integer of at least 2.")
}
if (is.na(roll_radius) || roll_radius < 0L) {
  stop("roll_radius must be a non-negative integer.")
}
if (!is.finite(tukey_alpha) || tukey_alpha <= 0 || tukey_alpha >= 1) {
  stop("tukey_alpha must be between 0 and 1.")
}
if (is.na(min_region_snps) || min_region_snps < 1L) {
  stop("min_region_snps must be a positive integer.")
}
if (!requireNamespace("bigutilsr", quietly = TRUE)) {
  stop(
    "R package 'bigutilsr' is required for the paper-matched ",
    "dist_ogk(), rollmean(), and tukey_mc_up() implementation."
  )
}
if (!requireNamespace("Cairo", quietly = TRUE)) {
  stop("R package 'Cairo' is required to write PCA LD-structure plots.")
}
if (!requireNamespace("ggplot2", quietly = TRUE) ||
    !requireNamespace("hexbin", quietly = TRUE)) {
  stop(
    "R packages 'ggplot2' and 'hexbin' are required for the ",
    "paper-style PCA loading plot."
  )
}

# These two functions intentionally keep the statistical core as close as
# possible to bigsnpr::bed_autoSVD(). Do not replace dist_ogk() with a classical
# covariance estimate or rollmean() with an unweighted moving average: either
# change would implement a different method from paper4-bedpca.
compute_robust_loading_distance <- function(pc_loadings) {
  # bigutilsr::dist_ogk() returns squared robust Mahalanobis distances.
  # bed_autoSVD() takes the square root before smoothing.
  sqrt(bigutilsr::dist_ogk(pc_loadings))
}

gaussian_smooth_by_chromosome <- function(
    statistic,
    chromosome,
    radius = 50L) {
  chromosome_indices <- split(seq_along(statistic), chromosome)
  smoothed <- rep(NA_real_, length(statistic))

  for (chromosome_name in names(chromosome_indices)) {
    index <- chromosome_indices[[chromosome_name]]
    if ((2L * radius + 1L) > length(index)) {
      stop(
        "Chromosome ", chromosome_name, " has only ", length(index),
        " PCA variants; Gaussian radius ", radius,
        " requires at least ", 2L * radius + 1L, " variants."
      )
    }
    smoothed[index] <- bigutilsr::rollmean(statistic[index], radius)
  }

  smoothed
}

output_dir <- dirname(output_prefix)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

message("Reading PCA SNP loadings from: ", loadings_file)
loadings <- utils::read.delim(
  loadings_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

metadata_columns <- c(
  "CHR", "POS", "SNP", "LOADING_ALLELE", "OTHER_ALLELE"
)
report_pc_columns <- paste0("PC", seq_len(n_report_pcs))
plot_pc_columns <- paste0("PC", seq_len(20L))
required_columns <- unique(c(
  metadata_columns,
  report_pc_columns,
  plot_pc_columns
))
missing_columns <- setdiff(required_columns, colnames(loadings))
if (length(missing_columns) > 0L) {
  stop(
    "PCA SNP-loading file is missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}
if (nrow(loadings) == 0L) {
  stop("PCA SNP-loading file contains no variants.")
}
if (anyNA(loadings[, metadata_columns, drop = FALSE])) {
  stop("PCA SNP-loading metadata contains missing values.")
}
if (anyDuplicated(loadings$SNP)) {
  stop("PCA SNP-loading file contains duplicate SNP IDs.")
}

loading_matrix <- as.matrix(
  loadings[, report_pc_columns, drop = FALSE]
)
storage.mode(loading_matrix) <- "double"
if (any(!is.finite(loading_matrix))) {
  stop("PC1-PC", n_report_pcs, " loadings contain non-finite values.")
}
if (nrow(loading_matrix) <= ncol(loading_matrix)) {
  stop("Too few variants for robust Mahalanobis distance estimation.")
}

plot_loading_matrix <- as.matrix(
  loadings[, plot_pc_columns, drop = FALSE]
)
storage.mode(plot_loading_matrix) <- "double"
if (any(!is.finite(plot_loading_matrix))) {
  stop("PC1-PC20 plot loadings contain non-finite values.")
}

# PCA loadings generated in 01b are autosomal. Accept either "1" or "chr1"
# spelling, then enforce chromosomes 1-22 so ordering and smoothing are stable.
chr_text <- sub("^CHR", "", toupper(trimws(as.character(loadings$CHR))))
chr_numeric <- suppressWarnings(as.integer(chr_text))
if (anyNA(chr_numeric) || any(!chr_numeric %in% seq_len(22L))) {
  stop("CHR must contain autosomes 1-22, optionally prefixed with 'chr'.")
}

position <- suppressWarnings(as.numeric(loadings$POS))
if (any(!is.finite(position)) || any(position < 1)) {
  stop("POS must contain positive finite genomic positions.")
}

# bed_autoSVD assumes variants are in genomic order before chromosome-wise
# smoothing. Sort explicitly so all 01b PCA branches behave identically.
variant_order <- order(chr_numeric, position, as.character(loadings$SNP))
loadings <- loadings[variant_order, , drop = FALSE]
loading_matrix <- loading_matrix[variant_order, , drop = FALSE]
plot_loading_matrix <- plot_loading_matrix[variant_order, , drop = FALSE]
chr_numeric <- chr_numeric[variant_order]
position <- position[variant_order]

message(
  "Computing robust Mahalanobis statistic from ",
  paste(report_pc_columns, collapse = ", "), "."
)
# dist_ogk() returns squared robust Mahalanobis distances. bed_autoSVD()
# applies sqrt() before Gaussian smoothing, so reproduce that exactly.
robust_mahalanobis <- compute_robust_loading_distance(loading_matrix)
if (length(robust_mahalanobis) != nrow(loadings) ||
    any(!is.finite(robust_mahalanobis))) {
  stop("Robust Mahalanobis calculation returned invalid values.")
}

message(
  "Applying chromosome-wise Gaussian smoothing with radius ",
  roll_radius, " variants."
)
chromosome_indices <- split(seq_len(nrow(loadings)), chr_numeric)
smoothed_statistic <- gaussian_smooth_by_chromosome(
  robust_mahalanobis,
  chr_numeric,
  roll_radius
)
if (any(!is.finite(smoothed_statistic))) {
  stop("Gaussian smoothing returned non-finite values.")
}

# With coef=NULL, tukey_mc_up() corrects for medcouple-estimated skewness and
# multiple testing. This matches bed_autoSVD() rather than standard Tukey 1.5.
threshold <- as.numeric(bigutilsr::tukey_mc_up(
  smoothed_statistic,
  alpha = tukey_alpha
))
if (length(threshold) != 1L || !is.finite(threshold)) {
  stop("Modified Tukey threshold is not a single finite value.")
}
is_outlier <- smoothed_statistic > threshold

# Regroup consecutive outliers within each chromosome. A candidate region must
# contain at least min_region_snps consecutive PCA variants, matching the
# bed_autoSVD() int.min.size rule.
region_id <- rep(NA_character_, nrow(loadings))
region_rows <- list()
region_counter <- 0L
for (chromosome in names(chromosome_indices)) {
  chromosome_index <- chromosome_indices[[chromosome]]
  chromosome_outliers <- chromosome_index[is_outlier[chromosome_index]]
  if (length(chromosome_outliers) == 0L) {
    next
  }

  run_group <- cumsum(c(TRUE, diff(chromosome_outliers) != 1L))
  outlier_runs <- split(chromosome_outliers, run_group)
  for (run in outlier_runs) {
    if (length(run) < min_region_snps) {
      next
    }
    region_counter <- region_counter + 1L
    current_region_id <- sprintf("LRLD%03d", region_counter)
    region_id[run] <- current_region_id
    peak_index <- run[which.max(smoothed_statistic[run])]
    region_rows[[region_counter]] <- data.frame(
      REGION_ID = current_region_id,
      CHR = as.integer(chromosome),
      START = min(position[run]),
      END = max(position[run]),
      N_CONSECUTIVE_OUTLIER_SNPS = length(run),
      FIRST_SNP = as.character(loadings$SNP[run[1L]]),
      LAST_SNP = as.character(loadings$SNP[run[length(run)]]),
      PEAK_SNP = as.character(loadings$SNP[peak_index]),
      PEAK_SMOOTHED_STAT = smoothed_statistic[peak_index],
      THRESHOLD = threshold,
      stringsAsFactors = FALSE
    )
  }
}

if (length(region_rows) > 0L) {
  regions <- do.call(rbind, region_rows)
} else {
  regions <- data.frame(
    REGION_ID = character(),
    CHR = integer(),
    START = numeric(),
    END = numeric(),
    N_CONSECUTIVE_OUTLIER_SNPS = integer(),
    FIRST_SNP = character(),
    LAST_SNP = character(),
    PEAK_SNP = character(),
    PEAK_SMOOTHED_STAT = numeric(),
    THRESHOLD = numeric(),
    stringsAsFactors = FALSE
  )
}

statistics <- data.frame(
  CHR = chr_numeric,
  POS = position,
  SNP = as.character(loadings$SNP),
  LOADING_ALLELE = as.character(loadings$LOADING_ALLELE),
  OTHER_ALLELE = as.character(loadings$OTHER_ALLELE),
  loadings[, report_pc_columns, drop = FALSE],
  ROBUST_MAHALANOBIS = robust_mahalanobis,
  GAUSSIAN_SMOOTHED_STAT = smoothed_statistic,
  TUKEY_THRESHOLD = rep(threshold, nrow(loadings)),
  IS_OUTLIER = is_outlier,
  REGION_ID = region_id,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

write_gzip_table <- function(x, filename) {
  connection <- gzfile(filename, open = "wt")
  on.exit(close(connection))
  utils::write.table(
    x,
    file = connection,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
}

statistics_file <- paste0(output_prefix, ".snp_stats.tsv.gz")
outlier_file <- paste0(output_prefix, ".outlier_snps.tsv.gz")
regions_file <- paste0(output_prefix, ".regions.tsv")
summary_file <- paste0(output_prefix, ".summary.tsv")
note_file <- paste0(output_prefix, ".note.txt")
statistic_plot_file <- paste0(output_prefix, ".smoothed_statistic.png")
loadings_plot_file <- paste0(output_prefix, ".PC1-PC20_loadings.png")

write_gzip_table(statistics, statistics_file)
write_gzip_table(statistics[is_outlier, , drop = FALSE], outlier_file)
utils::write.table(
  regions,
  file = regions_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

summary_table <- data.frame(
  METRIC = c(
    "INPUT_LOADING_FILE",
    "PCS_USED",
    "N_VARIANTS",
    "ROLL_RADIUS_VARIANTS",
    "TUKEY_ALPHA",
    "TUKEY_THRESHOLD",
    "N_OUTLIER_VARIANTS",
    "MIN_REGION_SNPS",
    "N_CANDIDATE_LRLD_REGIONS",
    "BIGUTILSR_VERSION"
  ),
  VALUE = c(
    normalizePath(loadings_file),
    paste(report_pc_columns, collapse = ","),
    nrow(statistics),
    roll_radius,
    format(tukey_alpha, scientific = FALSE),
    format(threshold, digits = 16),
    sum(is_outlier),
    min_region_snps,
    nrow(regions),
    as.character(utils::packageVersion("bigutilsr"))
  ),
  stringsAsFactors = FALSE
)
utils::write.table(
  summary_table,
  file = summary_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

note_lines <- c(
  "DEEP PCA long-range LD structure diagnostic",
  "",
  "IMPORTANT:",
  "This is a QC-only, single-round report. It does not remove variants,",
  "recompute PCA, or modify pca.eigenvec, pca.eigenval, genotype data,",
  "or pca.loadings.tsv.gz.",
  "",
  "Method matched to the outlier-detection step in bigsnpr::bed_autoSVD():",
  paste0(
    "1. Use SNP loadings from ",
    paste(report_pc_columns, collapse = ", "), "."
  ),
  "2. Compute sqrt(bigutilsr::dist_ogk(loadings)).",
  paste0(
    "3. Apply bigutilsr::rollmean separately within each chromosome ",
    "using a radius of ", roll_radius, " variants (a ",
    2L * roll_radius + 1L, "-variant Gaussian-weighted window)."
  ),
  paste0(
    "4. Compute bigutilsr::tukey_mc_up(smoothed statistic, alpha = ",
    tukey_alpha, ")."
  ),
  "   With coef=NULL, this adjusts for medcouple skewness and multiple testing.",
  paste0(
    "5. Flag variants above the threshold and report runs of at least ",
    min_region_snps, " consecutive flagged variants as candidate regions."
  ),
  "",
  paste0("Modified Tukey threshold: ", format(threshold, digits = 16)),
  paste0("Outlier variants: ", sum(is_outlier)),
  paste0("Candidate long-range LD regions: ", nrow(regions)),
  "The modified Tukey threshold is one global value calculated from all",
  "chromosome-wise smoothed SNP statistics; it is not SNP-specific.",
  "",
  "Interpretation:",
  paste0(
    "The robust Mahalanobis distance combines each SNP's loadings across ",
    paste(report_pc_columns, collapse = ", "),
    " into one statistic. It is large when a SNP has an unusual joint loading",
    " pattern after robustly accounting for the centre, scale, and covariance",
    " of the loading columns."
  ),
  paste0(
    "Gaussian smoothing replaces each SNP's raw distance with a weighted local",
    " average over up to ", 2L * roll_radius + 1L,
    " variants on the same chromosome. Variants nearest the focal SNP receive",
    " the largest weights. This emphasises clusters of unusual SNPs and",
    " suppresses isolated spikes."
  ),
  "IS_OUTLIER marks every SNP whose smoothed statistic exceeds the modified",
  paste0(
    "Tukey threshold. REGION_ID is assigned only to runs of at least ",
    min_region_snps, " consecutive outlier SNPs."
  ),
  "A candidate region may represent long-range LD structure, but its reported",
  "start and end are smoothing-based QC boundaries, not fine-mapped biological",
  "boundaries.",
  "The PC1-PC20 loading figure follows the paper's Figure 1 plotting method:",
  "SNP column index versus loading is displayed with ggplot2::geom_hex(), and",
  "the viridis fill colour represents the number of SNPs in each hexagonal bin.",
  "Grey vertical lines mark chromosome boundaries; x-axis labels are",
  "chromosome numbers placed at each chromosome midpoint.",
  "It is not automatically removed because this 01b QC step is report-only.",
  "Review the tables and plots before deciding whether PCA should be recomputed",
  "after excluding candidate long-range LD variants."
)
writeLines(note_lines, con = note_file)

# The LD-statistic plot uses genomic variant index with chromosome labels.
genome_index <- seq_len(nrow(statistics))
chromosome_starts <- vapply(chromosome_indices, min, integer(1))
chromosome_ends <- vapply(chromosome_indices, max, integer(1))
chromosome_midpoints <- (chromosome_starts + chromosome_ends) / 2
chromosome_label_order <- as.character(seq_len(22L))
chromosome_axis_breaks <- chromosome_midpoints[chromosome_label_order]
chromosome_axis_breaks <- chromosome_axis_breaks[
  !is.na(chromosome_axis_breaks)
]
chromosome_axis_labels <- names(chromosome_axis_breaks)

# Always include the single global threshold in the y-axis. Without this,
# datasets with zero outliers have threshold > max(smoothed_statistic), which
# places the threshold line outside the default plotting range.
statistic_y_range <- range(c(smoothed_statistic, threshold), finite = TRUE)
statistic_y_padding <- diff(statistic_y_range) * 0.05
if (!is.finite(statistic_y_padding) || statistic_y_padding == 0) {
  statistic_y_padding <- max(1, abs(statistic_y_range)) * 0.05
}
statistic_y_limits <- statistic_y_range +
  c(-statistic_y_padding, statistic_y_padding)
outlier_count <- sum(is_outlier)

Cairo::CairoPNG(
  filename = statistic_plot_file,
  width = 2600,
  height = 1100,
  pointsize = 20,
  res = 200
)
graphics::par(mar = c(7, 5, 3, 1))
graphics::plot(
  genome_index,
  smoothed_statistic,
  pch = 16,
  cex = 0.22,
  col = "#3B6FB6",
  xaxt = "n",
  xlab = "",
  ylab = "Smoothed robust distance",
  ylim = statistic_y_limits,
  main = paste0(
    "PCA loading LD-structure diagnostic (",
    paste0("PC1-PC", n_report_pcs), ")"
  ),
  sub = paste0(
    "One global modified Tukey threshold = ",
    format(threshold, digits = 6),
    "; SNPs above threshold = ", outlier_count
  )
)
graphics::axis(
  1,
  at = chromosome_axis_breaks,
  labels = FALSE,
  tick = FALSE
)
plot_user_limits <- graphics::par("usr")
graphics::text(
  x = chromosome_axis_breaks,
  y = plot_user_limits[3] - diff(plot_user_limits[3:4]) * 0.035,
  labels = chromosome_axis_labels,
  srt = 45,
  adj = c(1, 1),
  xpd = NA
)
graphics::mtext("Chromosome", side = 1, line = 5.4)
graphics::abline(v = chromosome_ends[-length(chromosome_ends)],
                 col = "grey85", lwd = 0.8)
graphics::abline(h = threshold, col = "#D73027", lwd = 2.5, lty = 2)
if (outlier_count > 0L) {
  graphics::points(
    genome_index[is_outlier],
    smoothed_statistic[is_outlier],
    pch = 16,
    cex = 0.32,
    col = "#D73027"
  )
  graphics::legend(
    "topright",
    legend = c(
      "SNP",
      paste0("SNP above threshold (n=", outlier_count, ")"),
      "One global modified Tukey threshold"
    ),
    col = c("#3B6FB6", "#D73027", "#D73027"),
    pch = c(16, 16, NA),
    lty = c(NA, NA, 2),
    bty = "n",
    cex = 0.8
  )
} else {
  graphics::legend(
    "topright",
    legend = c(
      "SNP",
      "One global modified Tukey threshold (no SNPs above)"
    ),
    col = c("#3B6FB6", "#D73027"),
    pch = c(16, NA),
    lty = c(NA, 2),
    bty = "n",
    cex = 0.8
  )
}
grDevices::dev.off()

# Reproduce the plotting approach used for Figure 1 in paper4-bedpca:
# one loading panel per PC, geom_hex() over SNP column index, a viridis count
# scale, and five columns of panels. The paper plotted PC1-PC40; DEEP plots
# PC1-PC20, producing four rows of five panels.
loading_plot_data <- data.frame(
  COLUMN_INDEX = rep(genome_index, times = length(plot_pc_columns)),
  LOADING = as.vector(
    plot_loading_matrix[, plot_pc_columns, drop = FALSE]
  ),
  PC = factor(
    rep(
      paste0("Loadings of ", plot_pc_columns),
      each = length(genome_index)
    ),
    levels = paste0("Loadings of ", plot_pc_columns)
  )
)
paper_style_loading_plot <- ggplot2::ggplot(
  loading_plot_data,
  ggplot2::aes(x = COLUMN_INDEX, y = LOADING)
) +
  ggplot2::geom_hex() +
  ggplot2::geom_vline(
    xintercept = chromosome_ends[-length(chromosome_ends)],
    colour = "grey70",
    linewidth = 0.3
  ) +
  ggplot2::facet_wrap(
    ggplot2::vars(PC),
    ncol = 5,
    scales = "free_y"
  ) +
  ggplot2::scale_fill_viridis_c(name = "SNP count") +
  ggplot2::scale_x_continuous(
    breaks = unname(chromosome_axis_breaks),
    labels = chromosome_axis_labels
  ) +
  ggplot2::labs(x = "Chromosome (SNP column index order)", y = NULL) +
  ggplot2::theme_minimal(base_size = 13) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1,
      vjust = 0.5,
      size = 8
    ),
    strip.text = ggplot2::element_text(face = "bold"),
    legend.position = "right"
  )

Cairo::CairoPNG(
  filename = loadings_plot_file,
  width = 5000,
  height = 3000,
  pointsize = 16,
  res = 200
)
print(paper_style_loading_plot)
grDevices::dev.off()

message(
  "PCA LD-structure QC complete: ", sum(is_outlier),
  " outlier variants and ", nrow(regions),
  " candidate long-range LD regions."
)
message("Report prefix: ", output_prefix)

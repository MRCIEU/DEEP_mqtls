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
# PCA, and repeats. DEEP runs this script in 01e as a QC report only, after 01c
# and 01d have already consumed the formal 01b PCs. It therefore MUST NOT alter
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
    if ((2L * radius + 1L) >= length(index)) {
      stop(
        "Chromosome ", chromosome_name, " has only ", length(index),
        " PCA variants; Gaussian radius ", radius,
        " requires more than ", 2L * radius + 1L, " variants."
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
required_columns <- c(metadata_columns, report_pc_columns)
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
loadings_plot_file <- paste0(output_prefix, ".PC1-PC10_loadings.png")

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
  "It is not automatically removed because 01e runs after downstream steps have",
  "already consumed the formal PCs from 01b. Review the tables and plots before",
  "deciding whether PCA should be recomputed in 01b."
)
writeLines(note_lines, con = note_file)

# Use genomic variant index, as in the paper's loading plots. Chromosome labels
# are placed at chromosome midpoints and boundaries are shown with grey lines.
genome_index <- seq_len(nrow(statistics))
chromosome_starts <- vapply(chromosome_indices, min, integer(1))
chromosome_ends <- vapply(chromosome_indices, max, integer(1))
chromosome_midpoints <- (chromosome_starts + chromosome_ends) / 2

Cairo::CairoPNG(
  filename = statistic_plot_file,
  width = 2600,
  height = 1100,
  pointsize = 20,
  res = 200
)
graphics::par(mar = c(5, 5, 3, 1))
graphics::plot(
  genome_index,
  smoothed_statistic,
  pch = 16,
  cex = 0.22,
  col = ifelse(is_outlier, "#D73027", "#3B6FB6"),
  xaxt = "n",
  xlab = "Chromosome",
  ylab = "Smoothed robust distance",
  main = paste0(
    "PCA loading LD-structure diagnostic (",
    paste0("PC1-PC", n_report_pcs), ")"
  )
)
graphics::axis(
  1,
  at = chromosome_midpoints,
  labels = names(chromosome_midpoints),
  tick = FALSE
)
graphics::abline(h = threshold, col = "#D73027", lwd = 2, lty = 2)
graphics::abline(v = chromosome_ends[-length(chromosome_ends)],
                 col = "grey85", lwd = 0.8)
graphics::legend(
  "topright",
  legend = c("Variant", "Above modified Tukey threshold", "Threshold"),
  col = c("#3B6FB6", "#D73027", "#D73027"),
  pch = c(16, 16, NA),
  lty = c(NA, NA, 2),
  bty = "n",
  cex = 0.8
)
grDevices::dev.off()

Cairo::CairoPNG(
  filename = loadings_plot_file,
  width = 2600,
  height = 3200,
  pointsize = 16,
  res = 200
)
graphics::par(
  mfrow = c(5, 2),
  mar = c(3.2, 4.2, 2.2, 0.8),
  oma = c(2, 1, 2, 0)
)
for (pc_column in report_pc_columns) {
  graphics::plot(
    genome_index,
    loading_matrix[, pc_column],
    pch = 16,
    cex = 0.18,
    col = ifelse(is_outlier, "#D73027", "#33333355"),
    xaxt = "n",
    xlab = "",
    ylab = "Loading",
    main = pc_column
  )
  graphics::axis(
    1,
    at = chromosome_midpoints,
    labels = names(chromosome_midpoints),
    tick = FALSE,
    cex.axis = 0.65
  )
  graphics::abline(h = 0, col = "grey70", lwd = 0.8)
  graphics::abline(v = chromosome_ends[-length(chromosome_ends)],
                   col = "grey90", lwd = 0.6)
}
graphics::mtext(
  "PCA SNP loadings; red variants exceed the smoothed robust-distance threshold",
  side = 3,
  outer = TRUE,
  line = 0.3,
  cex = 1.1
)
graphics::mtext("Chromosome", side = 1, outer = TRUE, line = 0.3)
grDevices::dev.off()

message(
  "PCA LD-structure QC complete: ", sum(is_outlier),
  " outlier variants and ", nrow(regions),
  " candidate long-range LD regions."
)
message("Report prefix: ", output_prefix)

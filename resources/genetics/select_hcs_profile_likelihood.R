#!/usr/bin/env Rscript

# Select a cohort-specific HC cutoff using the Zhu-Ghodsi profile-likelihood
# elbow criterion. The complete HC file is preserved; a separate file
# containing FID, IID and HC1..HCk is written for later downstream use.
#
# Args:
#   1. HC spectrum CSV (HC, index, prop_explained)
#   2. Complete HC file (FID IID HC1..HCn)
#   3. Selected HC output file
#   4. Per-cutoff profile-likelihood CSV
#   5. Selection summary CSV
#   6. Profile-likelihood PNG

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop(
    "Usage: select_hcs_profile_likelihood.R ",
    "<spectrum_csv> <complete_hcs_file> <selected_hcs_file> ",
    "<likelihood_csv> <selection_csv> <plot_png>"
  )
}

spectrum_file <- args[1]
complete_hcs_file <- args[2]
selected_hcs_file <- args[3]
likelihood_file <- args[4]
selection_file <- args[5]
plot_file <- args[6]

if (normalizePath(complete_hcs_file, mustWork = FALSE) ==
    normalizePath(selected_hcs_file, mustWork = FALSE)) {
  stop("The selected HC output must not overwrite the complete HC file")
}
if (file.exists(selected_hcs_file)) {
  file.remove(selected_hcs_file)
}

required_spectrum_columns <- c("HC", "index", "prop_explained")
minimum_noise_hcs <- 5L

profile_likelihood <- function(values, analysis, index_offset = 0L) {
  n_values <- length(values)
  if (n_values < 3L) {
    stop(
      analysis,
      " profile likelihood requires at least 3 component values; found ",
      n_values
    )
  }

  result <- rbindlist(lapply(seq_len(n_values - 1L), function(q) {
    signal <- values[seq_len(q)]
    noise <- values[(q + 1L):n_values]
    signal_mean <- mean(signal)
    noise_mean <- mean(noise)
    signal_sse <- sum((signal - signal_mean)^2)
    noise_sse <- sum((noise - noise_mean)^2)
    pooled_variance <- (signal_sse + noise_sse) / (n_values - 2L)

    if (!is.finite(pooled_variance) || pooled_variance <= 0) {
      log_likelihood <- NA_real_
    } else {
      log_likelihood <- sum(
        dnorm(
          signal,
          mean = signal_mean,
          sd = sqrt(pooled_variance),
          log = TRUE
        )
      ) + sum(
        dnorm(
          noise,
          mean = noise_mean,
          sd = sqrt(pooled_variance),
          log = TRUE
        )
      )
    }

    data.table(
      analysis = analysis,
      local_q = q,
      cutoff_hc = q + index_offset,
      n_signal = length(signal),
      n_noise = length(noise),
      signal_mean = signal_mean,
      noise_mean = noise_mean,
      signal_sse = signal_sse,
      noise_sse = noise_sse,
      pooled_variance = pooled_variance,
      log_likelihood = log_likelihood
    )
  }))

  if (!any(is.finite(result$log_likelihood))) {
    stop(
      analysis,
      " profile likelihood has no finite candidate values; ",
      "the pooled variance is zero or non-finite"
    )
  }

  selected_row <- which.max(result$log_likelihood)
  result[, selected := FALSE]
  result[selected_row, selected := TRUE]
  result
}

message("Reading HC spectrum: ", spectrum_file)
spectrum <- fread(spectrum_file)

missing_columns <- setdiff(required_spectrum_columns, names(spectrum))
if (length(missing_columns) > 0L) {
  stop(
    "HC spectrum is missing required column(s): ",
    paste(missing_columns, collapse = ", ")
  )
}

if (nrow(spectrum) < 4L) {
  stop(
    "At least 4 HCs are required to estimate both primary and secondary ",
    "profile-likelihood elbows; found ",
    nrow(spectrum)
  )
}

expected_index <- seq_len(nrow(spectrum))
expected_hc <- paste0("HC", expected_index)
if (!identical(as.integer(spectrum$index), expected_index) ||
    !identical(as.character(spectrum$HC), expected_hc)) {
  stop("HC spectrum must contain consecutive HC1..HCn rows in index order")
}

values <- as.numeric(spectrum$prop_explained)
if (anyNA(values) || any(!is.finite(values))) {
  stop("prop_explained contains NA or non-finite values")
}
if (any(values <= 0)) {
  stop("prop_explained values must all be greater than zero")
}
if (any(diff(values) > 0)) {
  stop("prop_explained values must be ordered from largest to smallest")
}

message("Calculating primary profile likelihood using HC1..HC", length(values))
primary <- profile_likelihood(
  values,
  analysis = "primary_full_spectrum",
  index_offset = 0L
)

message("Calculating secondary profile likelihood using HC2..HC", length(values))
secondary <- profile_likelihood(
  values[-1L],
  analysis = "secondary_excluding_HC1",
  index_offset = 1L
)

likelihood_results <- rbindlist(list(primary, secondary), use.names = TRUE)
fwrite(likelihood_results, likelihood_file)

primary_selected <- primary[selected == TRUE]
secondary_selected <- secondary[selected == TRUE]
final_cutoff <- as.integer(secondary_selected$cutoff_hc)
selection_status <- if (secondary_selected$n_noise < minimum_noise_hcs) {
  "review_required"
} else {
  "ok"
}

selection_summary <- data.table(
  n_hcs_available = length(values),
  primary_cutoff = as.integer(primary_selected$cutoff_hc),
  secondary_local_q = as.integer(secondary_selected$local_q),
  final_cutoff = final_cutoff,
  minimum_noise_hcs = minimum_noise_hcs,
  remaining_noise_hcs = as.integer(secondary_selected$n_noise),
  input_scale = "raw_prop_explained",
  status = selection_status
)
fwrite(selection_summary, selection_file)

plot_data <- copy(likelihood_results)
plot_data[, analysis_label := factor(
  analysis,
  levels = c("primary_full_spectrum", "secondary_excluding_HC1"),
  labels = c(
    "Primary profile likelihood: HC1 onwards",
    "Secondary profile likelihood: HC1 excluded"
  )
)]

selection_lines <- plot_data[selected == TRUE, .(
  analysis_label,
  cutoff_hc
)]

profile_plot <- ggplot(
  plot_data,
  aes(x = cutoff_hc, y = log_likelihood)
) +
  geom_line(linewidth = 0.5, na.rm = TRUE) +
  geom_point(size = 1, na.rm = TRUE) +
  geom_vline(
    data = selection_lines,
    aes(xintercept = cutoff_hc),
    colour = "red",
    linetype = "dashed",
    linewidth = 0.6
  ) +
  facet_wrap(~analysis_label, ncol = 1, scales = "free_y") +
  labs(
    title = "Profile-likelihood selection of Haplotype Components",
    subtitle = paste0(
      "Primary cutoff: HC",
      primary_selected$cutoff_hc,
      "; final secondary cutoff: HC",
      final_cutoff
    ),
    x = "Cutoff after HC",
    y = "Profile log-likelihood"
  ) +
  theme_bw()

ggsave(
  filename = plot_file,
  plot = profile_plot,
  width = 10,
  height = 8,
  dpi = 300
)

if (selection_status == "review_required") {
  message(
    "Profile-likelihood cutoff HC",
    final_cutoff,
    " leaves only ",
    secondary_selected$n_noise,
    " noise HCs; at least ",
    minimum_noise_hcs,
    " are required. Selection requires review, so no selected HC file ",
    "was written."
  )
} else {
  message("Reading complete HC file: ", complete_hcs_file)
  complete_hcs <- fread(complete_hcs_file)
  required_hc_columns <- paste0("HC", seq_len(final_cutoff))
  required_output_columns <- c("FID", "IID", required_hc_columns)
  missing_hc_columns <- setdiff(required_output_columns, names(complete_hcs))

  if (length(missing_hc_columns) > 0L) {
    stop(
      "Complete HC file is missing required column(s): ",
      paste(missing_hc_columns, collapse = ", ")
    )
  }

  selected_hcs <- complete_hcs[, ..required_output_columns]
  fwrite(
    selected_hcs,
    selected_hcs_file,
    sep = " ",
    col.names = TRUE
  )
  message(
    "Selected HC file written with HC1..HC",
    final_cutoff,
    ": ",
    selected_hcs_file
  )
}

message("Profile-likelihood results written to: ", likelihood_file)
message("Profile-likelihood selection written to: ", selection_file)
message("Profile-likelihood plot written to: ", plot_file)

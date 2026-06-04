#!/usr/bin/env Rscript

# Scatter plots of HC pairs (HC1 vs HC2, HC2 vs HC3, HC1 vs HC3) coloured by
# population / geographic factors from the winsorized phenotype file.
# The three pairs are drawn side by side on a single figure (one page per
# population factor, sharing a common legend).
# Intended as a QC check for population structure captured by HCs.
#
# Args:
#   1. hcs_file              space-separated HCs file (FID IID HC1..HCn)
#   2. covariates_file       combined covariates table (header; *_factor cols)
#   3. pheno_rdata           winsorized phenotype RData (loads object `pheno`)
#   4. study_name            used in output filenames
#   5. out_dir               output directory

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(viridis)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
    stop("Usage: hc_pop_scatter.R <hcs_file> <covariates_file> <winsorized_phenotype_file> <study_name> <out_dir>")
}
hcs_file        <- args[1]
covariates_file <- args[2]
pheno_rdata     <- args[3]
study_name      <- args[4]
out_dir         <- args[5]

candidate_factors <- c(
    "Birth_place_factor",
    "Population_group_factor",
    "Location_factor",
    "Urban_rural_factor",
    "Language_factor"
)

hcs <- fread(hcs_file)

covs <- as.data.table(
    read.table(covariates_file, header = TRUE, stringsAsFactors = FALSE,
               colClasses = c(Sex_factor = "character"))
)

# Any *_factor column is categorical by convention, even when its values are
# purely numeric (e.g. Plate_factor = 1, 2, 3). Force them to factor so they
# use the discrete colour scale. Keep only IID + the factor columns and
# Age_numeric from the covariates file (HCs provide the coordinates; cov is
# only used for colours).
cov_factor_cols <- grep("_factor$", names(covs), value = TRUE)
for (fc in cov_factor_cols) covs[[fc]] <- as.factor(covs[[fc]])
cov_keep <- c("IID", cov_factor_cols, intersect("Age_numeric", names(covs)))
covs <- covs[, cov_keep, with = FALSE]

env <- new.env()
load(pheno_rdata, envir = env)
if (!"pheno" %in% ls(env)) {
    stop("Expected object `pheno` in ", pheno_rdata, " (found: ",
         paste(ls(env), collapse = ", "), ")")
}
pheno <- as.data.table(env$pheno)

# From the phenotype data take only the population-structure factors, and
# only those not already supplied by the covariates file (avoids duplicates).
pheno_factors <- setdiff(intersect(candidate_factors, names(pheno)),
                         names(covs))
pheno <- pheno[, c("IID", pheno_factors), with = FALSE]
for (fc in pheno_factors) pheno[[fc]] <- as.factor(pheno[[fc]])

hcs[, IID := as.character(IID)]
covs[, IID := as.character(IID)]
pheno[, IID := as.character(IID)]

dat <- merge(merge(hcs, covs, by = "IID"), pheno, by = "IID")
cat("Merged samples for HC population scatter:", nrow(dat), "\n")
if (nrow(dat) == 0) {
    stop("No samples remain after merging HCs, covariates and phenotype data.")
}

# Colour variables: Age_numeric (continuous) plus every *_factor column
# available (all factors from the covariates file plus the population-structure
# factors from the phenotype).
present_colour_vars <- c(intersect("Age_numeric", names(dat)),
                         grep("_factor$", names(dat), value = TRUE))
if (length(present_colour_vars) == 0) {
    cat("No colour variables available; producing uncoloured scatter only.\n")
}

all_hc_pairs <- list(c("HC1", "HC2"), c("HC2", "HC3"), c("HC1", "HC3"))

# Keep only pairs where both HC columns are present in the merged data.
hc_pairs <- Filter(function(pr) all(pr %in% names(dat)), all_hc_pairs)
skipped_pairs <- Filter(function(pr) !all(pr %in% names(dat)), all_hc_pairs)
if (length(skipped_pairs) > 0) {
    cat("Skipping HC pair(s) (column not found):",
        paste(sapply(skipped_pairs, paste, collapse = "v"), collapse = ", "), "\n")
}
if (length(hc_pairs) == 0) {
    stop("No HC pairs available to plot.")
}

pair_levels <- sapply(hc_pairs, function(pr) paste0(pr[1], " vs ", pr[2]))

# Reshape data into long format so HC pairs can be faceted in one plot.
# Each pair becomes a data.table with identical column names (hc_x, hc_y,
# pair [, colour factor]) so binding is trivial.
reshape_pairs <- function(d, pairs, factor_cols) {
    rbindlist(lapply(pairs, function(pr) {
        out <- data.table(
            hc_x = d[[pr[1]]],
            hc_y = d[[pr[2]]],
            pair = factor(paste0(pr[1], " vs ", pr[2]), levels = pair_levels)
        )
        for (fc in factor_cols) out[[fc]] <- d[[fc]]
        out
    }))
}

out_pdf <- file.path(out_dir,
                     sprintf("HCs_population_structure_%s.pdf", study_name))

pdf(out_pdf, width = 15, height = 6)

if (length(present_colour_vars) == 0) {
    long <- reshape_pairs(dat, hc_pairs, character(0))
    p <- ggplot(long, aes(x = hc_x, y = hc_y)) +
        geom_point(alpha = 0.7, size = 1.4) +
        facet_wrap(~ pair, nrow = 1, scales = "free") +
        labs(title = "HC scatter (HC1 vs HC2, HC2 vs HC3, HC1 vs HC3)",
             subtitle = paste0("Study: ", study_name, "  |  N = ", nrow(dat),
                               "  |  no population factors available"),
             x = NULL, y = NULL) +
        theme_bw()
    print(p)
} else {
    for (cvar in present_colour_vars) {
        d <- dat[!is.na(get(cvar))]
        if (nrow(d) == 0) {
            cat("  ", cvar, ": no non-missing samples, skipped.\n", sep = "")
            next
        }
        long <- reshape_pairs(d, hc_pairs, cvar)
        if (is.numeric(long[[cvar]])) {
            p <- ggplot(long, aes(x = hc_x, y = hc_y,
                                  colour = .data[[cvar]])) +
                geom_point(alpha = 0.7, size = 1.4) +
                facet_wrap(~ pair, nrow = 1, scales = "free") +
                scale_colour_viridis() +
                labs(title = paste0("HC scatter coloured by ", cvar),
                     subtitle = paste0("Study: ", study_name, "  |  N = ",
                                       nrow(d), "  |  continuous"),
                     x = NULL, y = NULL, colour = cvar) +
                theme_bw()
            print(p)
            next
        }
        long[[cvar]] <- as.factor(long[[cvar]])
        nlev <- nlevels(long[[cvar]])
        p <- ggplot(long, aes(x = hc_x, y = hc_y, colour = .data[[cvar]])) +
            geom_point(alpha = 0.7, size = 1.4) +
            facet_wrap(~ pair, nrow = 1, scales = "free") +
            scale_colour_viridis(discrete = TRUE) +
            labs(title = paste0("HC scatter coloured by ", cvar),
                 subtitle = paste0("Study: ", study_name, "  |  N = ", nrow(d),
                                   "  |  ", nlev, " levels"),
                 x = NULL, y = NULL, colour = cvar) +
            theme_bw()
        if (nlev > 12) {
            # Too many levels for a usable legend: keep the colours but
            # drop the legend entirely.
            p <- p + guides(colour = "none")
        } else if (nlev > 8) {
            p <- p + theme(legend.text = element_text(size = 6),
                           legend.key.size = unit(0.3, "cm")) +
                guides(colour = guide_legend(
                    ncol = 2, override.aes = list(size = 2)))
        }
        print(p)
    }
}

dev.off()
cat("Wrote ", out_pdf, "\n", sep = "")

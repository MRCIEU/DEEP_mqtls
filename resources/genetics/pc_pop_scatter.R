#!/usr/bin/env Rscript

# Scatter plots of genetic PC pairs (PC1 vs PC2, PC2 vs PC3, PC1 vs PC3)
# coloured by Age_numeric, Sex_factor and population / geographic factors.
# Continuous variables (e.g. Age_numeric) use a continuous viridis scale;
# categorical variables use a discrete viridis scale.
# The three pairs are drawn side by side on a single figure (one page per
# colour variable). Intended as a QC check for population structure captured
# by the genetic PCs.
#
# Args:
#   1. covariates_file       combined covariates table (header; contains
#                            IID, genetic_pc1.., Age_numeric, Sex_factor, ...)
#   2. pheno_rdata           winsorized phenotype RData (loads object `pheno`)
#   3. study_name            used in output filenames
#   4. out_dir               output directory

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(viridis)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
    stop("Usage: pc_pop_scatter.R <covariates_file> <winsorized_phenotype_file> <study_name> <out_dir>")
}
covariates_file <- args[1]
pheno_rdata     <- args[2]
study_name      <- args[3]
out_dir         <- args[4]

covs <- as.data.table(
    read.table(covariates_file, header = TRUE, stringsAsFactors = FALSE,
               colClasses = c(Sex_factor = "character"))
)

# Any *_factor column is categorical by convention, even when its values are
# purely numeric (e.g. Plate_factor = 1, 2, 3). Force them to factor so they
# use the discrete colour scale rather than being treated as continuous.
cov_factor_cols <- grep("_factor$", names(covs), value = TRUE)
for (fc in cov_factor_cols) covs[[fc]] <- as.factor(covs[[fc]])

candidate_factors <- c(
    "Birth_place_factor",
    "Population_group_factor",
    "Location_factor",
    "Urban_rural_factor",
    "Language_factor"
)

env <- new.env()
load(pheno_rdata, envir = env)
if (!"pheno" %in% ls(env)) {
    stop("Expected object `pheno` in ", pheno_rdata, " (found: ",
         paste(ls(env), collapse = ", "), ")")
}
pheno <- as.data.table(env$pheno)

# Take only IID + population factors present in the phenotype data, so the
# merge does not duplicate columns shared with the covariates file.
pheno_factors <- intersect(candidate_factors, names(pheno))
pheno <- pheno[, c("IID", pheno_factors), with = FALSE]
for (fc in pheno_factors) pheno[[fc]] <- as.factor(pheno[[fc]])

covs[, IID := as.character(IID)]
pheno[, IID := as.character(IID)]

dat <- merge(covs, pheno, by = "IID")
cat("Merged samples for PC population scatter:", nrow(dat), "\n")
if (nrow(dat) == 0) {
    stop("No samples remain after merging covariates and phenotype data.")
}

# Colour variables: age, plus every *_factor column available (all population
# / technical factors from the covariates file and the phenotype data). Using
# the *_factor suffix means newly added factor columns are picked up
# automatically without editing this script.
candidate_colour_vars <- unique(c(
    "Age_numeric",
    grep("_factor$", names(dat), value = TRUE)
))
present_colour_vars <- intersect(candidate_colour_vars, names(dat))
missing_colour_vars <- setdiff(candidate_colour_vars, names(dat))
if (length(missing_colour_vars) > 0) {
    cat("Skipping colour variables (not present):",
        paste(missing_colour_vars, collapse = ", "), "\n")
}
if (length(present_colour_vars) == 0) {
    cat("No colour variables available; producing uncoloured scatter only.\n")
}

all_pc_pairs <- list(c("genetic_pc1", "genetic_pc2"),
                     c("genetic_pc2", "genetic_pc3"),
                     c("genetic_pc1", "genetic_pc3"))

# Keep only pairs where both PC columns are present in the merged data.
pc_pairs <- Filter(function(pr) all(pr %in% names(dat)), all_pc_pairs)
skipped_pairs <- Filter(function(pr) !all(pr %in% names(dat)), all_pc_pairs)
if (length(skipped_pairs) > 0) {
    cat("Skipping PC pair(s) (column not found):",
        paste(sapply(skipped_pairs, paste, collapse = "v"), collapse = ", "), "\n")
}
if (length(pc_pairs) == 0) {
    stop("No PC pairs available to plot.")
}

# Pretty pair labels (PC1 vs PC2 rather than genetic_pc1 vs genetic_pc2).
pair_label <- function(pr) {
    paste0(sub("genetic_pc", "PC", pr[1]), " vs ",
           sub("genetic_pc", "PC", pr[2]))
}
pair_levels <- sapply(pc_pairs, pair_label)

# Reshape data into long format so PC pairs can be faceted in one plot.
# Each pair becomes a data.table with identical column names (pc_x, pc_y,
# pair [, colour variable]). The colour column keeps its original type so
# numeric variables stay numeric (continuous colour scale).
reshape_pairs <- function(d, pairs, colour_col) {
    rbindlist(lapply(pairs, function(pr) {
        out <- data.table(
            pc_x = d[[pr[1]]],
            pc_y = d[[pr[2]]],
            pair = factor(pair_label(pr), levels = pair_levels)
        )
        if (length(colour_col) > 0) out[[colour_col]] <- d[[colour_col]]
        out
    }))
}

out_pdf <- file.path(out_dir,
                     sprintf("PCs_population_structure_%s.pdf", study_name))

pdf(out_pdf, width = 15, height = 6)

if (length(present_colour_vars) == 0) {
    long <- reshape_pairs(dat, pc_pairs, character(0))
    p <- ggplot(long, aes(x = pc_x, y = pc_y)) +
        geom_point(alpha = 0.7, size = 1.4) +
        facet_wrap(~ pair, nrow = 1, scales = "free") +
        labs(title = "Genetic PC scatter (PC1 vs PC2, PC2 vs PC3, PC1 vs PC3)",
             subtitle = paste0("Study: ", study_name, "  |  N = ", nrow(dat),
                               "  |  no colour variables available"),
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
        long <- reshape_pairs(d, pc_pairs, cvar)
        is_num <- is.numeric(long[[cvar]])
        p <- ggplot(long, aes(x = pc_x, y = pc_y, colour = .data[[cvar]])) +
            geom_point(alpha = 0.7, size = 1.4) +
            facet_wrap(~ pair, nrow = 1, scales = "free") +
            theme_bw()
        if (is_num) {
            nval <- length(unique(long[[cvar]]))
            p <- p + scale_colour_viridis() +
                labs(title = paste0("Genetic PC scatter coloured by ", cvar),
                     subtitle = paste0("Study: ", study_name, "  |  N = ",
                                       nrow(d), "  |  continuous"),
                     x = NULL, y = NULL, colour = cvar)
        } else {
            long[[cvar]] <- as.factor(long[[cvar]])
            nlev <- nlevels(long[[cvar]])
            p <- ggplot(long, aes(x = pc_x, y = pc_y,
                                  colour = .data[[cvar]])) +
                geom_point(alpha = 0.7, size = 1.4) +
                facet_wrap(~ pair, nrow = 1, scales = "free") +
                scale_colour_viridis(discrete = TRUE) +
                labs(title = paste0("Genetic PC scatter coloured by ", cvar),
                     subtitle = paste0("Study: ", study_name, "  |  N = ",
                                       nrow(d), "  |  ", nlev, " levels"),
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
        }
        print(p)
    }
}

dev.off()
cat("Wrote ", out_pdf, "\n", sep = "")

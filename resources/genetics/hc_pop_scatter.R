#!/usr/bin/env Rscript

# Scatter plots of HC pairs (HC1 vs HC2, HC3 vs HC4) coloured by
# population / geographic factors from the winsorized phenotype file.
# Intended as a QC check for population structure captured by HCs.
#
# Args:
#   1. hcs_file              space-separated HCs file (FID IID HC1..HCn)
#   2. pheno_rdata           winsorized phenotype RData (loads object `pheno`)
#   3. study_name            used in output filenames
#   4. out_dir               output directory

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
    stop("Usage: hc_pop_scatter.R <hcs_file> <winsorized_phenotype_file> <study_name> <out_dir>")
}
hcs_file    <- args[1]
pheno_rdata <- args[2]
study_name  <- args[3]
out_dir     <- args[4]

hcs <- fread(hcs_file)

env <- new.env()
load(pheno_rdata, envir = env)
if (!"pheno" %in% ls(env)) {
    stop("Expected object `pheno` in ", pheno_rdata, " (found: ",
         paste(ls(env), collapse = ", "), ")")
}
pheno <- as.data.table(env$pheno)

hcs[, IID := as.character(IID)]
pheno[, IID := as.character(IID)]

dat <- merge(hcs, pheno, by = "IID")
cat("Merged samples for HC population scatter:", nrow(dat), "\n")
if (nrow(dat) == 0) {
    stop("No samples remain after merging HCs and phenotype data.")
}

candidate_factors <- c(
    "Birth_place_factor",
    "Population_group_factor",
    "Location_factor",
    "Urban_rural_factor",
    "Language_factor"
)
present_factors <- intersect(candidate_factors, names(dat))
missing_factors <- setdiff(candidate_factors, names(dat))
if (length(missing_factors) > 0) {
    cat("Skipping (not present in phenotype data):",
        paste(missing_factors, collapse = ", "), "\n")
}
if (length(present_factors) == 0) {
    cat("No population/geographic factors available; producing uncoloured scatter only.\n")
}

plot_hc_pair <- function(xvar, yvar, out_pdf) {
    if (!all(c(xvar, yvar) %in% names(dat))) {
        cat("Skipping ", xvar, " vs ", yvar,
            ": column(s) missing from HCs file.\n", sep = "")
        return(invisible(NULL))
    }
    pdf(out_pdf, width = 7, height = 6)
    if (length(present_factors) == 0) {
        p <- ggplot(dat, aes_string(x = xvar, y = yvar)) +
            geom_point(alpha = 0.7, size = 1.4) +
            labs(title = paste0(xvar, " vs ", yvar),
                 subtitle = paste0("Study: ", study_name,
                                   "  |  N = ", nrow(dat),
                                   "  |  no population factors available"),
                 x = xvar, y = yvar) +
            theme_bw()
        print(p)
        dev.off()
        cat("Wrote ", out_pdf, "\n", sep = "")
        return(invisible(NULL))
    }
    for (fac in present_factors) {
        d <- dat[!is.na(get(fac))]
        if (nrow(d) == 0) {
            cat("  ", fac, ": no non-missing samples, skipped.\n", sep = "")
            next
        }
        d[[fac]] <- as.factor(d[[fac]])
        nlev <- nlevels(d[[fac]])
        p <- ggplot(d, aes_string(x = xvar, y = yvar, colour = fac)) +
            geom_point(alpha = 0.7, size = 1.4) +
            labs(title = paste0(xvar, " vs ", yvar, " coloured by ", fac),
                 subtitle = paste0("Study: ", study_name,
                                   "  |  N = ", nrow(d),
                                   "  |  ", nlev, " levels"),
                 x = xvar, y = yvar, colour = fac) +
            theme_bw()
        if (nlev > 12) {
            p <- p + theme(legend.position = "none") +
                labs(subtitle = paste0("Study: ", study_name,
                                       "  |  N = ", nrow(d),
                                       "  |  ", nlev,
                                       " levels (legend hidden)"))
        }
        print(p)
    }
    dev.off()
    cat("Wrote ", out_pdf, "\n", sep = "")
}

plot_hc_pair("HC1", "HC2",
             file.path(out_dir,
                       sprintf("HCs_population_structure_HC1HC2_%s.pdf",
                               study_name)))
plot_hc_pair("HC3", "HC4",
             file.path(out_dir,
                       sprintf("HCs_population_structure_HC3HC4_%s.pdf",
                               study_name)))

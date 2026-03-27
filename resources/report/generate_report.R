#!/usr/bin/env Rscript

# Generate summary QC report for DEEP pipeline 01 modules
# This script renders the Rmd template to HTML
suppressPackageStartupMessages({
    library(rmarkdown)
    library(knitr)
    library(kableExtra)
    library(ggplot2)
    library(magick)
    library(base64enc)
})

args <- commandArgs(TRUE)
home_directory     <- as.character(args[1])
scripts_directory  <- as.character(args[2])
section_01_dir     <- as.character(args[3])
study_name         <- as.character(args[4])
ancestry           <- as.character(args[5])
genome_build       <- as.character(args[6])
meth_chip          <- as.character(args[7])
n_pcs              <- as.numeric(args[8])
n_hcs              <- as.numeric(args[9])

# Define paths
output_dir <- file.path(section_01_dir, "report")
rmd_file <- file.path(scripts_directory, "resources", "report", "summary_report.Rmd")
output_file <- file.path(output_dir, "summary_report.html")

# Check Rmd exists
if (!file.exists(rmd_file)) {
    stop("Rmd template not found: ", rmd_file)
}
# Create output directory if needed
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}
# Render the report
message("Rendering report...")
message("  Input:  ", rmd_file)
message("  Output: ", output_file)
rmarkdown::render(
    input = rmd_file,
    output_file = basename(output_file),
    output_dir = output_dir,
    params = list(
        home_directory = home_directory,
        scripts_directory = scripts_directory,
        section_01_dir = section_01_dir,
        study_name = study_name,
        ancestry = ancestry,
        genome_build = genome_build,
        meth_chip = meth_chip,
        n_pcs = n_pcs,
        n_hcs = n_hcs
    ),
    quiet = FALSE
)
message("Report generated: ", output_file)
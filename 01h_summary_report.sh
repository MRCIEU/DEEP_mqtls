#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

exec &> >(tee ${section_01h_logfile})

print_version
echo "==============================================="
echo "01h: Generating Summary Report"
echo "==============================================="
echo ""

# Check that required previous steps have been completed
echo "Checking previous step outputs..."
# Required files from 01a
if [ ! -f "${cohort_descriptives}" ]; then
    echo "ERROR: cohort_descriptives.RData not found. Please run 01a first."
    exit 1
fi

# Required files from 01b
if [ ! -f "${pcaplot}" ]; then
    echo "WARNING: pcaplot.pdf not found from 01b."
fi

# Required files from 01c
if [ ! -f "${cohort_descriptives_was}" ]; then
    echo "WARNING: cohort_descriptives_was.RData not found from 01c."
fi

# Create output directory
mkdir -p "${section_01_dir}/report"
echo ""
echo "Generating HTML summary report..."
echo ""
${R_directory}Rscript resources/report/generate_report.R \
    "${home_directory}" \
    "${scripts_directory}" \
    "${section_01_dir}" \
    "${study_name}" \
    "${ancestry}" \
    "${genome_build}" \
    "${meth_chip}" \
    "${n_pcs}" \
    "${n_hcs}"

report_file="${section_01_dir}/report/summary_report.html"

if [ -f "${report_file}" ]; then
    echo ""
    echo "==============================================="
    echo "Summary report generated successfully!"
    echo "Location: ${report_file}"
    echo "==============================================="
else
    echo "ERROR: Report generation failed."
    exit 1
fi

echo ""
echo "Successfully completed script 01h"
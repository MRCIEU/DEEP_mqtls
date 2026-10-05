#!/bin/bash
# Generate Phase I QC HTML reports for extracted cohort results.
#
# Usage: bash 03_summary_report.sh <extracted_dir> <report_dir> [cohort ...]
#
#   extracted_dir  directory holding one subdirectory per cohort, each with
#                  a config file and results/01/
#   report_dir     where the HTML reports are written (created if missing)
#   cohort ...     optional: zero, one or more cohort names, separated by
#                  spaces. Each must match a subdirectory name in
#                  extracted_dir exactly (case-sensitive). If omitted, every
#                  subdirectory of extracted_dir is used. A cohort with no
#                  config file or results/01/ is skipped with a warning.
#
# Examples:
#   # All cohorts found in extracted_dir
#   bash 03_summary_report.sh sftp_downloads/extracted reports
#
#   # One cohort
#   bash 03_summary_report.sh sftp_downloads/extracted reports cohort_A
#
#   # Several cohorts
#   bash 03_summary_report.sh sftp_downloads/extracted reports cohort_A cohort_B

if [ $# -lt 2 ]; then
    sed -n '4,/^$/p' "$0" | sed 's/^# \{0,1\}//'
    exit 1
fi

EXTRACTED_DIR="${1%/}"
REPORT_DIR="${2%/}"
shift 2

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT="${SCRIPT_DIR}/generate_deep_qc_report_v2.py"

if [ ! -d "$EXTRACTED_DIR" ]; then
    echo "ERROR: extracted_dir not found: $EXTRACTED_DIR"
    exit 1
fi

if [ ! -f "$SCRIPT" ]; then
    echo "ERROR: report generator not found: $SCRIPT"
    exit 1
fi

mkdir -p "$REPORT_DIR"

if [ $# -gt 0 ]; then
    COHORTS=("$@")
else
    COHORTS=()
    for d in "$EXTRACTED_DIR"/*/; do
        [ -d "$d" ] && COHORTS+=("$(basename "$d")")
    done
fi

if [ ${#COHORTS[@]} -eq 0 ]; then
    echo "ERROR: no cohorts found in $EXTRACTED_DIR"
    exit 1
fi

for COHORT in "${COHORTS[@]}"; do
    CONFIG_FILE="${EXTRACTED_DIR}/${COHORT}/config"
    RESULTS_DIR="${EXTRACTED_DIR}/${COHORT}/results/01"
    REPORT_OUT="${REPORT_DIR}/${COHORT}_PhaseI_QC_report.html"

    echo "=== Generating report for $COHORT ==="

    if [ ! -f "$CONFIG_FILE" ]; then
        echo "WARNING: config file not found for $COHORT — skipping"
        continue
    fi

    if [ ! -d "$RESULTS_DIR" ]; then
        echo "WARNING: results/01 dir not found for $COHORT — skipping"
        continue
    fi

    if python3 "$SCRIPT" "$CONFIG_FILE" "$REPORT_OUT" --results "$RESULTS_DIR"; then
        echo "$COHORT report done → $REPORT_OUT"
    else
        echo "ERROR: Report generation failed for $COHORT"
    fi
done

echo "All reports done. Check $REPORT_DIR"

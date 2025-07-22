#!/bin/bash

# Check scripts up to date
git pull

# Initialize variables
config_file="./config"

# Activate the environment
if [ -z "$R_directory" ] && [ -z "$Python_directory" ] && [ -z "$Python2_directory" ]; then
    # Users are using system default R & Python, activate conda environment
    if ! command -v mamba &> /dev/null; then
        echo "mamba not found. Please install mamba first."
        exit 1
    fi
    # Initialize mamba shell for current session
    eval "$(mamba shell hook --shell bash)"
    mamba activate deep_mqtl
    echo "Current conda environment: $CONDA_DEFAULT_ENV"
else
    # Users have specified custom R/Python directories
    echo "Custom R/Python directories specified"
fi

# Parse options using getopts
while getopts "c:" opt; do
    case $opt in
        c) config_file=$OPTARG ;;
        *) echo "Usage: $0 -c <config_file>"
           exit 1 ;;
    esac
done

# Shift option arguments, so $1 becomes the first positional argument
shift $((OPTIND - 1))

# Initialize an empty string
concatenated=""

# Loop through all arguments
for arg in "$@"; do
    concatenated="$concatenated$arg "
done


set -e
echo "-----------------------------------------------"
echo ""
echo "Using config located at:" ${config_file}
echo ""
echo "-----------------------------------------------"

source ${config_file}
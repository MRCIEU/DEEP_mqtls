#!/bin/bash

# Check scripts up to date
git pull

# Initialize variables
config_file="./config"

# Activate the environment
if [ -z "$R_directory" ] && [ -z "$Python_directory" ] && [ -z "$Python2_directory" ]; then
    if command -v mamba &> /dev/null; then
        # Check if mamba is available and contains deep_env environment
        if mamba env list | awk 'NF > 0 && $1 !~ /^#/ && $1 !~ /^\// {print $1}' | grep -Fxq 'deep_env'; then
            echo "found deep_env environment in mamba"
            echo "Using mamba to run the script"
            RUN_CMD="mamba run -n deep_env"
        fi
    fi

    if [ -z "$RUN_CMD" ] && command -v conda &> /dev/null; then
        # Check if conda is available and contains deep_env environment
        if conda env list | awk 'NF > 0 && $1 !~ /^#/ && $1 !~ /^\// {print $1}' | grep -Fxq 'deep_env'; then
            echo "found deep_env environment in conda"
            echo "Using conda to run the script"
            RUN_CMD="conda run -n deep_env"
        fi
    fi

    if [ -z "$RUN_CMD" ]; then
        echo "ERROR: deep_env environment not found in mamba or conda."
        exit 1
    fi

else
    echo "Using specified Python3_directory"
    RUN_CMD=""
fi

export RUN_CMD

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
# generate pcs on the probes aross all the control probles
# quantile normalize the datasets together

# https://github.com/perishky/meffil/wiki/Functional-normalizing-separate-datasets

#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

exec &> >(tee ${section_01f_logfile})
print_version

${R_directory}Rscript resources/methylation/Func_normalization_bw_datasets.R \
		${study_name} \
		${normalized_between_datasets}

echo "Successfully completed script 01f"
# generate pcs on the probes aross all the control probles
# quantile normalize the datasets together

# https://github.com/perishky/meffil/wiki/Functional-normalizing-separate-datasets

#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

exec &> >(tee ${section_01f_logfile})
print_version

containsElement () {
	local e
	for e in "${@:2}"; do [[ "$e" == "$1" ]] && return 0; done
	echo "There is no method for ${1}."
	echo "Please run:"
	echo "./01f-normalization.sh [arg]"
	echo "where arg is an optional argument that can be one of:"
	printf '%s\n' ${@:2}
	return 1
}

arg=""
declare -a sections=('shrink' 'expand')


if [ -z "${1}" ]; then
    echo "Usage: ./01f-normalization.sh [arg]"
    echo "where arg is an optional argument that can be one of:"
    printf '%s\n' "${sections[@]}"
    exit 0
else
    arg="${1}"
    containsElement "${1}" "${sections[@]}"
fi

section_message () {

	echo "-----------------------------------------------"
	echo ""
	echo "$1 section"
	echo ""
	echo "to run this part on its own type:"
	echo "$ ./01f-normalization.sh $1"
	echo ""
	echo "-----------------------------------------------"
	echo ""
	echo ""

}

if [ "$arg" = "shrink" ]
then
	section_message "shrink"

	${R_directory}Rscript resources/methylation/Func_normalization_bw_datasets.R \
		0 \
		${meth_qc_obj} \
		${study_name} \
		${section_01_dir} \
		"shrink"

fi

if [ "$arg" = "expand" ]
then
	section_message "expand"

	${R_directory}Rscript resources/methylation/Func_normalization_bw_datasets.R \
		${betas} \
		${meth_qc_obj} \
		${study_name} \
		${home_directory} \
		"expand"
fi
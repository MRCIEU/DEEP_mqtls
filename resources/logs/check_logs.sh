#!/usr/bin/env bash
vercomp () {
	if [[ $1 == $2 ]]
	then
		echo "Correct script version"
		return 0
	fi
	local IFS=.
	local i ver1=($1) ver2=($2)
	# fill empty fields in ver1 with zeros
	for ((i=${#ver1[@]}; i<${#ver2[@]}; i++))
	do
		ver1[i]=0
	done
	for ((i=0; i<${#ver1[@]}; i++))
	do
		if [[ -z ${ver2[i]} ]]
		then
			# fill empty fields in ver2 with zeros
			ver2[i]=0
		fi
		if ((10#${ver1[i]} > 10#${ver2[i]}))
		then
			echo "Script version greater than required"
			return 0
		fi
		if ((10#${ver1[i]} < 10#${ver2[i]}))
		then
			echo ""
			echo "PROBLEM"
			echo "This analysis was performed on an outdated script."
			echo "Expecting at least version $2"
			echo "But the logs show that this was run on version $1"
			echo "Please run 'git pull' and then re-run the analysis."
			echo ""
			return 1
		fi
	done
	echo "Correct script version"
	return 0
}

compare_version () {

    section=$1
    suffix=${section: -1}
    log_dir="${section_01_dir}/logs_${suffix}"
    log_file=$(ls -1t "${log_dir}/log"* 2>/dev/null | head -n 1)

    if [ -z "$log_file" ]; then
        echo ""
        echo "WARNING"
        echo "No logfile found in ${log_dir}"
        echo "Please run 'git pull' and check that no updates were made to the ${section} script you are checking."
        echo "If updates were made then please re-run this ${section} script."
        echo ""
        return 0
    fi

    version_used=$(grep "DEEP version" "$log_file" | head -n 1 | sed 's/.*DEEP version \([0-9.]*\).*/\1/')
    if [ -z "$version_used" ]; then
        echo ""
        echo "WARNING"
        echo "No version number found in logfile: $log_file"
        echo "The scripts you used could be out of date."
        echo "Please run 'git pull' and check that no updates were made to the ${section} script you are checking."
        echo "If updates were made then please re-run this ${section} script."
        echo ""
        return 0
    fi

    version_required=$(grep "section_${section}" resources/logs/versions.txt | cut -d " " -f 2)
    echo "Version required: ${version_required}"
    echo "Version used:    ${version_used}"
    vercomp ${version_used} ${version_required}
}

check_latest_chunk_logs () {
    local section="$1" script_name="$2"
    shift 2
    local chunks=("$@")
    local log_dir="${section_01_dir}/logs_${section: -1}"
    local log_files=()
    local chunk log_file status selected_log legacy
    local success_count=0

    for log_file in "${log_dir}"/log*.txt; do
        if [ -f "$log_file" ]; then
            log_files+=("$log_file")
        fi
    done
    if [ ${#log_files[@]} -gt 0 ]; then
        # Keep the existing convention: most recently modified log first.
        mapfile -t log_files < <(LC_ALL=C ls -1t -- "${log_files[@]}")
    fi

    for chunk in "${chunks[@]}"; do
        selected_log=""
        status="absent"
        for log_file in "${log_files[@]}"; do
            legacy=0
            if [ "$section" = "01g" ] && [ "${log_file##*/}" = "log.txt" ]; then
                legacy=1
            fi
            if ! status=$(awk -v section="$section" -v chunk="$chunk" \
                -v legacy="$legacy" '
                BEGIN {
                    state = "absent"
                    start = chunk " section"
                    success = "Successfully completed script " section " " chunk " chunk"
                }
                { sub(/\r$/, "") }
                $0 == start {
                    state = "incomplete"
                    active = 1
                    next
                }
                /^[[:alnum:]_]+ section$/ { active = 0 }
                active && $0 == success { state = "complete" }
                legacy && $0 == "Successfully completed script 01g" { whole_success = 1 }
                END {
                    if (state == "absent" && whole_success)
                        print "legacy_complete"
                    else
                        print state
                }
            ' "$log_file"); then
                status="unreadable"
            fi
            if [ "$status" != "absent" ]; then
                selected_log="$log_file"
                break
            fi
        done

        if [ "$status" = "complete" ] || [ "$status" = "legacy_complete" ]; then
            success_count=$((success_count + 1))
            echo "${section} ${chunk}: completed successfully; logfile: ${selected_log}"
            if [ "$status" = "legacy_complete" ]; then
                echo "  Using the explicit whole-module completion marker in the legacy log."
            fi
        elif [ "$status" = "incomplete" ]; then
            echo "${section} ${chunk}: incomplete; logfile: ${selected_log}"
            echo "  Latest attempt started without its completion marker; older successes are ignored."
        elif [ "$status" = "unreadable" ]; then
            echo "${section} ${chunk}: cannot check logfile: ${selected_log}"
        else
            echo "${section} ${chunk}: no run record found in ${log_dir}"
        fi
    done

    echo "Successful chunks: $success_count/${#chunks[@]}"
    if [ $success_count -eq ${#chunks[@]} ]; then
        echo "${script_name} completed successfully."
    else
        echo "Problem: ${script_name} did not complete successfully ($success_count/${#chunks[@]} chunks)"
        exit 1
    fi
}

check_logs_01c () {
    check_latest_chunk_logs "01c" "01c-check_phenotypes_and_methylation.sh" \
        methy_outlier check_phenotype predict_age_smoking cell_counts ewas meth_pcs combine_covariates
}

check_logs_01g () {
    check_latest_chunk_logs "01g" "01g-HCs.sh" vcf hc gwas
}

check_logs_01 () {

	exec &> >(tee ${section_01a_uploadlog})

	compare_version "01a"
	log_files=("${section_01_dir}/logs_a/log"*)
    if [ ${#log_files[@]} -gt 0 ] && grep -i -q "success" "${log_files[@]}"; then
		echo "01a-check_data.sh completed successfully."
	else
		echo "Problem: 01a-check_data.sh did not complete successfully"
		exit 1
	fi

	compare_version "01b"
	if grep -i -q "success" ${section_01b_logfile}; then
		echo "01b-process_genetic_data.sh completed successfully."
	else
		echo "Problem: 01b-process_genetic_data.sh did not complete successfully"
		exit 1
	fi

	compare_version "01c"
	check_logs_01c
	
	compare_version "01d"
	if grep -i -q "success" ${section_01d_logfile}; then
		echo "01d-mqtl_controls.sh completed successfully."
	else
		echo "Problem: 01d-mqtl_controls.sh did not complete successfully"
		exit 1
	fi

	echo "Skip 01e log file check: 01e-genetic_pc_gwas.sh is no longer required."
	# compare_version "01e"
	# if grep -i -q "success" ${section_01e_logfile}; then
	# 	echo "01e-genetic_pc_gwas.sh completed successfully."
	# else
	# 	echo "Problem: 01e-genetic_pc_gwas.sh did not complete successfully"
	# 	exit 1
	# fi

	if [ -z "${idat_directory}" ]; then
    	echo "You don't have idat dir, skip 01f log file check. Please ensure you don't have access to idat files."
	else
		echo "IDAT directory found, checking 01f log file."
		compare_version "01f"
		if ls "${section_01_dir}/logs_f"/log*.txt >/dev/null 2>&1 && \
   			grep -i -q "Shrunk QC objects created" "${section_01_dir}/logs_f"/log*.txt; then
    		echo "01f-normalization.sh completed successfully."
		else
			echo "Problem: 01f-normalization.sh did not complete successfully"
			exit 1
		fi
	fi

	if [ -z "${vcf_dir}" ]; then
    	echo "You don't have vcf dir, skip 01g log file check. Please ensure you don't have phased genetic data."
	else
		echo "VCF directory found, checking 01g log file."
		compare_version "01g"
		check_logs_01g
	fi
}

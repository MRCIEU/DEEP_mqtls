#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

usage () {
	echo $"Usage: $0 [-c <config_file>] <pipeline section> {check|upload} [sftp|manual]"
}

checkFirstArg () {
	local e
	for e in "${@:2}"; do [[ "$e" == "$1" ]] && return 0; done
	echo $"Error: $1 is not a valid section identifier"
	echo $"Need to specify a value from 01"
	usage
	exit 1
}

checkSecondArg () {
	local e
	for e in "${@:2}"; do [[ "$e" == "$1" ]] && return 0; done
	echo $"Error: $2 is not a valid action"
	echo $"Specify either 'check' or 'upload'"
	usage
	exit 1
}

checkUploadMode () {
	local e
	for e in "${@:2}"; do [[ "$e" == "$1" ]] && return 0; done
	echo $"Error: $1 is not a valid upload mode"
	echo $"Specify either 'sftp' or 'manual'"
	usage
	exit 1
}

source resources/logs/check_logs.sh
source resources/logs/check_results.sh

sections=("01")
checkFirstArg "$1" "${sections[@]}"

actions=("check" "upload")
checkSecondArg "$2" "${actions[@]}"

upload_mode="${3:-sftp}"
upload_modes=("sftp" "manual")
if [[ "$2" = "upload" ]]; then
	checkUploadMode "$upload_mode" "${upload_modes[@]}"
elif [[ -n "$3" ]]; then
	echo $"Error: upload mode can only be specified with the 'upload' action"
	usage
	exit 1
fi

if [[ -n "$4" ]]; then
	echo $"Error: too many arguments"
	usage
	exit 1
fi

echo ""
echo "Checking log files for $1"
eval "check_logs_$1"

echo ""
echo "Checking results for $1"
eval "check_results_$1"

echo ""
echo "Section $1 has been successfully completed!"

port="-P 2222"

if [[ "$2" = "upload" && ( $1 = "01" ) ]]; then

	echo ""
	if [[ "$upload_mode" = "sftp" ]]; then
		if command -v sftp >/dev/null 2>&1; then
    echo "sftp detected"

    if ! sftp $port -oIdentityFile="$key" -oBatchMode=no -b - "${sftp_username}@${sftp_address}:${sftp_path}" <<EOF
bye
EOF
	then
        echo "sftp connection failed"
        exit 1
    fi

    echo "Connection established"
else
    echo "sftp is not available."
    echo "Please ensure you can connect to the server using sftp before tarring results."
	exit 1
fi
	else
		echo "Manual upload mode selected; skipping sftp connection check."
	fi


	echo ""
	echo "Tarring results, log and config files"
	mkdir -p ${home_directory}/results/config/
	if [[ "$config_file" = /* ]]; then
    	target="${scripts_directory}/$(basename "$config_file")"
    	if [[ "${config_file}" != "$target" ]]; then
        cp ${config_file} ${scripts_directory}/
    	fi
    	config_basename=$(basename "$config_file")
    	tar czf ${home_directory}/results/config/${study_name}_config.tar -C ${scripts_directory}/ ./${config_basename} ./resources/parameters
    	if [[ "${config_file}" != "$target" ]]; then
        	rm ${scripts_directory}/${config_basename}
    	fi
	else
    	tar czf ${home_directory}/results/config/${study_name}_config.tar -C ${scripts_directory}/ ./${config_file} ./resources/parameters 
	fi
    echo "Successfully created config archive"	

    echo "Generating md5 checksum"
    md5sum ${home_directory}/results/config/${study_name}_config.tar > ${home_directory}/results/config/${study_name}_config.md5sum
    echo "Encrypting files"
    gpg --output ${home_directory}/results/config/${study_name}_config.tar.aes --symmetric --cipher-algo AES256 ${home_directory}/results/config/${study_name}_config.tar
    echo ""
	
	suff="tgz"
	flags="czf"
    tar ${flags} ${home_directory}/results/${study_name}_${1}.${suff} -C ${home_directory} results/${1}
    echo "Generating md5 checksum"
    md5sum ${home_directory}/results/${study_name}_${1}.${suff} > ${home_directory}/results/${study_name}_${1}.md5sum
    echo "Encrypting files"
    gpg --output ${home_directory}/results/${study_name}_${1}.${suff}.aes --symmetric --cipher-algo AES256 ${home_directory}/results/${study_name}_${1}.${suff}
    echo ""

	if [[ "$upload_mode" = "sftp" ]]; then
		if command -v sftp >/dev/null 2>&1; then
        sftp $port -oIdentityFile=$key -oBatchMode=no -b - ${sftp_username}@${sftp_address}:${sftp_path} <<EOF
        dir
        cd ../upload
        put ${home_directory}/results/${study_name}_${1}.md5sum
        put ${home_directory}/results/${study_name}_$1.${suff}.aes
		put ${home_directory}/results/config/${study_name}_config.tar.aes
		put ${home_directory}/results/config/${study_name}_config.md5sum 
bye
EOF
		fi
	else
		echo "Manual upload mode selected; upload these files manually:"
		echo "${home_directory}/results/${study_name}_${1}.md5sum"
		echo "${home_directory}/results/${study_name}_${1}.${suff}.aes"
		echo "${home_directory}/results/config/${study_name}_config.tar.aes"
		echo "${home_directory}/results/config/${study_name}_config.md5sum"
	fi
fi

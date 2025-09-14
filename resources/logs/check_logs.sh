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

    logfile="section_${1}_logfile"
    version_used=$(grep "DEEP version" ${!logfile}* | head -n 1 | cut -d " " -f 3)
	if [ "${version_used}" = "" ]
	then
		echo ""
		echo "WARNING"
		echo "No version number found. You are probably running an old version of git."
		echo "The scripts you used could be out of date."
		echo "Please run 'git pull' and check that no updates were made to the ${1} script you are checking."
		echo "If updates were made then please re-run this ${1} script."
		echo ""
		return 0
	fi

	version_required=$(grep "section_${1}" resources/logs/versions.txt | cut -d " " -f 2)
	echo "Version required: ${version_required}"
	echo "Version used: ${version_used}"
	vercomp ${version_used} ${version_required}

}

check_logs_01 () {

	compare_version "01a"
	if grep -i -q "success" ${section_01_logfile}; then
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
	if grep -i -q "success" ${section_01c_logfile}; then
		echo "01c-check_phenotypes_and_methylation.sh completed successfully."
	else
		echo "Problem: 01c-check_phenotypes_and_methylation.sh did not complete successfully"
		exit 1
	fi

	compare_version "01d"
	if grep -i -q "success" ${section_01d_logfile}; then
		echo "01d-mqtl_controls.sh completed successfully."
	else
		echo "Problem: 01d-mqtl_controls.sh did not complete successfully"
		exit 1
	fi

	compare_version "01e"
	if grep -i -q "success" ${section_01e_logfile}; then
		echo "01e-genetic_pc_gwas.sh completed successfully."
	else
		echo "Problem: 01e-genetic_pc_gwas.sh did not complete successfully"
		exit 1
	fi
}
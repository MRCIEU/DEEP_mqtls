#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

checkFirstArg () {
	local e
	for e in "${@:2}"; do [[ "$e" == "$1" ]] && return 0; done
	echo $"Error: $1 is not a valid section identifier"
	echo $"Need to specify a value from 01"
	echo $"Usage: $0 <pipeline section> {check|upload}"
	exit 1
}

checkSecondArg () {
	local e
	for e in "${@:2}"; do [[ "$e" == "$1" ]] && return 0; done
	echo $"Error: $2 is not a valid action"
	echo $"Specify either 'check' or 'upload'"
	echo $"Usage: $0 <pipeline section> {check|upload}"
	exit 1
}

source resources/logs/check_logs.sh
source resources/logs/check_results.sh

sections=("01")
checkFirstArg "$1" "${sections[@]}"

actions=("check" "upload")
checkSecondArg "$2" "${actions[@]}"

echo ""
echo "Checking log files for $1"
eval "check_logs_$1"

echo ""
echo "Checking results for $1"
eval "check_results_$1"

echo ""
echo "Section $1 has been successfully completed!"

if [[ "$2" = "upload" && ( $1 = "01" ) ]]
then

	echo ""
	temp=`which sshpass 2>/dev/null | wc -l`
	port="-P 2222"
	if [[ ! "${temp}" = "0" &&  $1 != "08" ]]
	then
		echo "sshpass detected"
		sftp $port -oIdentityFile=$key -oBatchMode=no -b - ${sftp_username}@${sftp_address}:${sftp_path} << !
bye
!
	elif [[ ! "${temp}" = "0" && $1 = "08" ]]
	then
		read -s -p "Enter SFTP password: " mypassword
	    export SSHPASS=${mypassword}
		sshpass -e sftp -P 22 -oIdentityFile=$key -oBatchMode=no -b - ${sftp_username_inversions}@${sftp_address_inversions}:${sftp_path_inversions} << !
	
bye
!
		echo "Connection established"
	else
		echo "sshpass is not installed."
		echo "The results will now be archived, once that is done they will be uploaded to the server"
	fi

	echo ""
	echo "Tarring results, log and config files"
	mkdir -p ${home_directory}/results/config/
	if [[ "$config_file" = /* ]]; then
		cp ${config_file} ${scripts_directory}/
		config_basename=$(basename "$config_file")
		tar czf ${home_directory}/results/config/${study_name}_config.tar -C ${scripts_directory}/ ./${config_basename} ./resources/parameters
		rm ${scripts_directory}/${config_basename}
	else
		tar czf ${home_directory}/results/config/${study_name}_config.tar -C ${scripts_directory}/ ./${config_file} ./resources/parameters 
	fi	
    echo "Successfully created config archive"	

    echo "Generating md5 checksum"
    md5sum ${home_directory}/results/config/${study_name}_config.tar > ${home_directory}/results/config/${study_name}_config.md5sum
    echo "Encrypting files"
    gpg --output ${home_directory}/results/config/${study_name}_config.tar.aes --symmetric --cipher-algo AES256 ${home_directory}/results/config/${study_name}_config.tar
    echo ""

	if [ "${1}" = "17" ]
	then
		suff="tar"
		flags="cf"
	else
		suff="tgz"
		flags="czf"
	fi

    if [[ $1 = "03a" || $1 = "03d" ]]
    then
        tar ${flags} ${home_directory}/results/${study_name}_${1}.${suff} -C ${home_directory} results/03/
        echo "Successfully created results archives"
#        echo "Generating md5 checksum"
#        md5sum ${home_directory}/results/${study_name}_${1}.${suff} > ${home_directory}/results/${study_name}_${1}.md5sum
#        echo "Encrypting files"
#        gpg --output ${home_directory}/results/${study_name}_${1}.${suff}.aes --symmetric --cipher-algo AES256 ${home_directory}//results/${study_name}_${1}.${suff}
        echo ""
    elif [[ $1 = "07" || $1 = "04" ]]
    then
        echo "Tarring results have been generated"
    else
	    tar ${flags} ${home_directory}/results/${study_name}_${1}.${suff} -C ${home_directory} results/${1}
        echo "Successfully created results archives"
    fi
    
    echo "Generating md5 checksum"
    md5sum ${home_directory}/results/${study_name}_${1}.${suff} > ${home_directory}/results/${study_name}_${1}.md5sum
    echo "Encrypting files"
    gpg --output ${home_directory}/results/${study_name}_${1}.${suff}.aes --symmetric --cipher-algo AES256 ${home_directory}//results/${study_name}_${1}.${suff}
    echo ""


if [[ ! "${temp}" = "0"  ]]
then
echo "Detecting sshpass"

	if [[  $1 = "08" ]]
	then
		echo "Inversions"
		export SSHPASS=${mypassword}
	    sshpass -e sftp -P 22 -oIdentityFile=$key -oBatchMode=no -b - ${sftp_username_inversions}@${sftp_address_inversions}:${sftp_path_inversions} << !
   		dir
		cd ${sftp_username_inversions}
   		put ${home_directory}/results/${study_name}_${1}.md5sum
   		put ${home_directory}/results/${study_name}_$1.${suff}.aes
		put ${home_directory}/results/config/${study_name}_config.tar.aes
		put ${home_directory}/results/config/${study_name}_config.md5sum 
   		bye
!
    else
        sftp $port -oIdentityFile=$key -oBatchMode=no -b - ${sftp_username}@${sftp_address}:${sftp_path} << !
        dir
        cd ../upload
        put ${home_directory}/results/${study_name}_${1}.md5sum
        put ${home_directory}/results/${study_name}_$1.${suff}.aes
		put ${home_directory}/results/config/${study_name}_config.tar.aes
		put ${home_directory}/results/config/${study_name}_config.md5sum 
        bye
!
	fi
else

read -s -p "Ready to upload? Press enter to continue: " anykey

sftp $port -oIdentityFile=$key ${sftp_username}@${sftp_address}:${sftp_path} <<EOF

    dir
    cd ../upload
    put ${home_directory}/results/${study_name}_${1}.md5sum
    put ${home_directory}/results/${study_name}_$1.${suff}.aes
	put ${home_directory}/results/config/${study_name}_config.tar.aes
	put ${home_directory}/results/config/${study_name}_config.md5sum 
EOF

fi

fi
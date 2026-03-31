#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

port="-P 2222"

# Check sftp availability
if ! command -v sftp >/dev/null 2>&1; then
    echo "sftp is not available."
    echo "Please ensure you can connect to the server using sftp before tarring results."
    exit 1
fi

echo "sftp detected"

# Test connection
if ! sftp $port -oIdentityFile="$key" -oBatchMode=no -b - "${sftp_username}@${sftp_address}:${sftp_path}" <<EOF
bye
EOF
then
    echo "sftp connection failed"
    exit 1
fi
echo "Connection established"

# Tar config files
echo ""
echo "Tarring config files"
mkdir -p ${home_directory}/results/config/

if [[ "$config_file" = /* ]]; then
    tar czf ${home_directory}/results/config/${study_name}_config.tar \
        ${config_file} \
        ${scripts_directory}/resources/parameters
else
    tar czf ${home_directory}/results/config/${study_name}_config.tar \
        -C ${scripts_directory}/ ./${config_file} ./resources/parameters
fi
echo "Successfully created config archive"

echo "Generating md5 checksum"
md5sum ${home_directory}/results/config/${study_name}_config.tar > ${home_directory}/results/config/${study_name}_config.md5sum
echo "Encrypting files"
gpg --output ${home_directory}/results/config/${study_name}_config.tar.aes --symmetric --cipher-algo AES256 ${home_directory}/results/config/${study_name}_config.tar
echo ""

# Tar results for section 01
echo "Tarring results for 01"
suff="tgz"
tar czf ${home_directory}/results/${study_name}_01.${suff} -C ${home_directory} results/01
echo "Generating md5 checksum"
md5sum ${home_directory}/results/${study_name}_01.${suff} > ${home_directory}/results/${study_name}_01.md5sum
echo "Encrypting files"
gpg --output ${home_directory}/results/${study_name}_01.${suff}.aes --symmetric --cipher-algo AES256 ${home_directory}/results/${study_name}_01.${suff}
echo ""

# Upload everything
echo "Uploading files"
sftp $port -oIdentityFile=$key -oBatchMode=no -b - ${sftp_username}@${sftp_address}:${sftp_path} <<EOF
cd ../upload_debug
put ${home_directory}/results/${study_name}_01.md5sum
put ${home_directory}/results/${study_name}_01.${suff}.aes
put ${home_directory}/results/config/${study_name}_config.tar.aes
put ${home_directory}/results/config/${study_name}_config.md5sum
bye
EOF

echo ""
echo "Upload complete!"
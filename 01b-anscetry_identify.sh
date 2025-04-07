# assigning ancestry labels to participants

source resources/setup.sh "$@"
set -- $concatenated

exec &> >(tee ${section_01b_logfile})
print_version

